

from gcsnap.rich_console import RichConsole 
from gcsnap.configuration import Configuration 

def main():

    console = RichConsole()
    console.print_title()

    # 1. Parse configuration and arguments
    config = Configuration()
    config.parse_arguments()

    # 2. parse targets
    console.print_step('Parsing targets...')


    
if __name__ == '__main__':
    main()