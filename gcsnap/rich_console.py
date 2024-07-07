import argparse
# pip install rich
from rich.console import Console
from rich.text import Text
from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn
from contextlib import contextmanager

import logging
logger = logging.getLogger(__name__) # inherits configuration from main logger

class RichConsole():    
    def __init__(self):
        self.console = Console()

        # colors to use
        self.color_grey = 'gray78'
        self.color_blue = 'steel_blue1'
        self.color_purple = 'plum1'
        self.color_gold = 'light_goldenrod2'
        self.color_red = 'indian_red1'
        self.color_green = 'pale_green1'

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['console']  # Remove the console as it's not picklable
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.console = Console()  # Recreate the Console object        

    def print_title(self) -> None:
        color = self.color_blue
        self.print_line(color)
        self.console.print(Text('GCsnap', style=f'bold {color}'))
        self.console.print(Text('GCsnap is a python-based, local tool that generates ' + 
                        'interactive snapshots of conserved protein-coding genomic contexts.',
                        style=color))
        self.console.print(Text('Thanks for using it!  💙', style=color))        
        self.print_line(color)  
        logger.info('GCsnap started')

    def print_final(self) -> None:
        color = self.color_blue
        self.print_line(color)
        self.console.print(Text('🏁  GCsnap finished successfully! 🥳🎆💫', style=color))
        self.print_line(color)  
        logger.info('GCsnap finished successfully')                            

    def print_line(self, color: str) -> None:
        self.console.print(Text('---------------------------------------------', style=color))

    def print_error(self, message: str) -> None:
        self.console.print(Text('❌  Error {}.'.format(message), style=self.color_red))   
        logger.error(f'{message}')

    def print_warning(self, message: str) -> None:
        self.console.print(Text('  ⚠️  {} Check gcsnap.log'.format(message), style=self.color_gold))   
        logger.warning(f'{message}')        

    def print_skipped_step(self, message: str) -> None:
        self.console.print(Text(message, style=self.color_grey))
        logger.info(f'{message}')

    def print_info(self, message: str) -> None:
        self.console.print(Text('\t{}'.format(message), style=self.color_grey))  
        logger.info(f'{message}')        
         
    def print_hint(self, message: str) -> None:
        self.console.print(Text('  💡  {}'.format(message), style=self.color_gold))   
        logger.warning(f'{message}')          

    def print_done(self, message: str) -> None:
        # print done message, make first character of message lowercase            
        self.console.print(Text('✅  Done {}'.format(message[0].lower() + 
                                                      message[1:]), style=self.color_green))     
        logger.info(f'Done {message}')        

    @contextmanager
    def status(self, message: str):
        logger.info(f'{message}')  
        with self.console.status(Text('{} ...'.format(message), style=self.color_grey)):
            yield  
        self.print_done(message)
        
    @contextmanager
    def progress(self, message: str, total: int):
        logger.info(f'{message}')  
        with Progress(
            TextColumn("{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeRemainingColumn(),
            console=self.console
        ) as progress:
            task_id = progress.add_task('[{}]\t{}'.format(self.color_grey, message), total=total)
            yield progress, task_id        
        self.print_done(message)

    def print_help(self, parser) -> None:
        usage = 'GCsnap --targets <targets> [Optional arguments]'
        epilog = 'Example: GCsnap --targets PHOL_ECOLI A0A0U4VKN7_9PSED A0'

        usage = Text(parser.format_usage(), style=self.color_grey)
        epilog = Text('\n' + parser.epilog, style=self.color_grey)

        # Print usage
        self.console.print(usage)
        # print config.yaml message
        self.console.print(Text('If CLI arguments not sepcified, default value from config.yaml', style=self.color_green))
        self.console.print(Text('--overwrite-config to replace values with CLI Arguments', style=self.color_green))
        self.console.print(Text('Or change them manually in config.yaml\n', style=self.color_green))

        # Print --targets separately
        for action in parser._actions:
            if '--targets' in action.option_strings:
                self.print_argument(action)
                break

        # Optional arguments
        self.console.print(Text('\nOptional arguments:', style=self.color_grey))
        # Print arguments
        for action in parser._actions:
            # skip --help and --targets
            if '--targets' in action.option_strings or '--help' in action.option_strings:
                continue
            self.print_argument(action)

        # Print epilog
        self.console.print(epilog)
        
    def print_argument(self, action) -> None:
        option_strings = ", ".join(action.option_strings)
        metavar = action.metavar or ""
        argument_help = Text(f"{option_strings} {metavar}", style=self.color_blue)
        if action.help:
            argument_help.append(f"\n    {action.help}", style=self.color_grey)

        if action.choices:
            choices_str = ', '.join(map(str, action.choices))
            argument_help.append(f'\n    Supported values: [{choices_str}]', style=self.color_grey)        

        if action.default != argparse.SUPPRESS:
            default_value = f"\n    [config.yaml value: {action.default}]"
            argument_help.append(default_value, style=self.color_green)

        self.console.print(argument_help)

