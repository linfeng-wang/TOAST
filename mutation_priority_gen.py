import csv
import argparse

def generate_csv_from_positions(position_string, output_file):
    """
    Generates a CSV file based on a list of comma-separated positions provided by the user.

    Parameters:
    - position_string (str): Comma-separated string of positions (e.g., "322168,553767,1077188").
    - output_file (str): Path to the output CSV file.
    """
    # Split the string into a list of positions
    positions = [int(pos.strip()) for pos in position_string.split(",")]

    # Create a list of dictionaries with genome positions
    data = []
    for i, pos in enumerate(positions, start=1):
        data.append({
            "sample_id": f"sample_{i}",
            "genome_pos": pos,
            "gene": f"gene_{i}",
            "change": f"change_{i}",
            "freq": 1,
            "type": "-",
            "sublin": "-",
            "drtype": "-",
            "drugs": "-",
            "weight": 1
        })

    # Write the data to a CSV file
    with open(output_file, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=data[0].keys())
        writer.writeheader()
        writer.writerows(data)

    print(f"CSV file generated: {output_file}")


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate a CSV file from a list of genome positions.")
    parser.add_argument(
        "-p", "--positions",
        required=True,
        help="Comma-separated list of genome positions (e.g., '322168,553767,1077188')."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file name (e.g., 'output.csv')."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Generate the CSV file
    generate_csv_from_positions(args.positions, args.output)


if __name__ == "__main__":
    main()


# python mutation_priority_gen.py -p '12204,12555,811492,811642,811753,1135798,1136017,1227217,1227518,4238120,4238963' -o /mnt/storage10/lwang/Projects/TOAST/cache/lin_snp01-29.csv