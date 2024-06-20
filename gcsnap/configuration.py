import os
import argparse

# pip install rich
from rich.console import Console
# pip install pyyaml
import yaml
# pip install ruamel.yaml
from ruamel.yaml.main import round_trip_dump as yaml_dump

from gcsnap.rich_console import RichConsole 


class CustomArgumentParser(argparse.ArgumentParser):
    def __init__(self, usage=None, epilog=None):
        super().__init__(usage=usage, epilog=epilog, add_help=False)

    def add_argument_from_config(self, config):
        for key, details in config.items():
            value_type = eval(details['type'])
            help_text = details['help']
            default_value = details['value']
            
            if value_type == bool:
                # argparse does not support bool type in config, so we need to convert it
                self.add_argument(f'--{key}', type=self.str_to_bool, default=default_value, help=help_text)
            else:
                self.add_argument(f'--{key}', type=value_type, default=default_value, help=help_text)

    @staticmethod
    def str_to_bool(value):
        if isinstance(value, bool):
            return value
        if value.lower() in ('yes', 'true', 'True' , 't', 'y', '1'):
            return True
        elif value.lower() in ('no', 'false', 'False', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    @staticmethod
    def str_to_type(type_str):
        if type_str == 'bool':
            return bool
        elif type_str == 'int':
            return int
        elif type_str == 'float':
            return float
        elif type_str == 'str':
            return str
        else:
            raise ValueError(f"Unsupported type: {type_str}")                


class Configuration:    
    def __init__(self):                     
        self.path = os.path.dirname(__file__)
        self.console = RichConsole()

        self.read_configuration_yaml()
        self.create_argument_parser()

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['parser']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.create_argument_parser()
        
    def read_configuration_yaml(self) -> None:   
        with open(os.path.join(self.path,'config.yaml'), 'r') as file:
            self.arguments = yaml.load(file, Loader=yaml.FullLoader)
            
    def write_configuration_yaml(self) -> None:
        out = self.arguments
        header_comment = (
            'Configuration file to set arguments for GCsnap.\n'
            'To change argument, change: value: entry.\n'
            'E.g. value: 1 to value: 2\n'
            '---------------------------------------\n')
        
        with open(os.path.join(self.path,'config.yaml'), 'w') as file:
            file.write('# ' + header_comment.replace('\n', '\n# ') + '\n')
            # 4 is the indentation space
            yaml_dump(self.arguments, default_flow_style=False, indent=4, block_seq_indent=4)
        
    def parse_arguments(self) -> None:
        args = self.parser.parse_args()
        if hasattr(args, 'help') and args.help:
            self.console.print_help(self.parser)
            exit(0)

        # Ensure 'targets' argument is present
        if args.targets is None:
            self.console.print_error('The following argument is required: --targets')
            exit(0) 

        # Update self.arguments dictionary with parsed values
        for arg in vars(args):
            if arg in self.arguments:
                self.arguments[arg]['value'] = getattr(args, arg)
            else:
                self.targets = getattr(args, 'targets')

        # handle arguments that require special treatment
        self.handle_special_arguments()

        # Write updated configuration to file
        if self.arguments['overwrite_config']['value']:
            self.write_configuration_yaml()
        
    def create_argument_parser(self) -> None:
        usage = 'GCsnap --targets <targets> [Optional arguments]'
        epilog = 'Example: GCsnap --targets PHOL_ECOLI A0A0U4VKN7_9PSED A0'
        
        self.parser = CustomArgumentParser(usage=usage, epilog=epilog)

        # help argument, as we overwrote the default --help with rich printing
        self.parser.add_argument('-h', '--help', action='store_true', 
                                 help='Show this help message and exit')
        # consequences, we can't use required=True but check it in parse_arguments()

        # required arguments
        self.parser.add_argument('-t', '--targets', nargs='+', type=str, 
                                 help='List of input targets. Can be a list of fasta files,' + 
                                 'a list of text files encompassing a list of protein sequence identifiers,' + 
                                 'a list of protein sequence identifiers, or a mix of them')
        
        # Add arguments from YAML
        self.parser.add_argument_from_config(self.arguments)

    def handle_special_arguments(self) -> None:
        # handle special arguments that require additional processing

        # TODO: Check how this with clans file works, seems to be a special case
        
        # clanse file
        targets = self.targets
        clans_file = self.arguments['clans_file']['value']
        if len(targets) == 1 and targets[0].endswith('.clans') and clans_file is None:
            clans_file = targets[0]
        if clans_file is not None:
            clans_file = os.path.abspath(clans_file)

        self.arguments['clans_file']['value'] = clans_file