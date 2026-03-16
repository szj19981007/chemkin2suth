import argparse
import sys
from .converter import convert_transport_file

def main():
    parser = argparse.ArgumentParser(
        prog="chemkin2suth",
        description="Convert Chemkin transport files to OpenFOAM Sutherland format."
    )
    parser.add_argument(
        "trans_input_file",
        help="Chemkin transport input file"    
    )
    parser.add_argument(
        "thermo_input_file",
        help="Chemkin thermo input file"    
    )
    parser.add_argument(
        "output_file",
        help="Output OpenFOAM transport dictionary"
    )

    args = parser.parse_args()

    try:
        convert_transport_file(
            args.trans_input_file,
            args.thermo_input_file,
            args.output_file
        )
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)