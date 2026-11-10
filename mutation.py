import random
import statistics
from Bio.Seq import Seq
import matplotlib.pyplot as plt

AMINO_ACIDS = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]


def mutate_position(sequence, position):
    # Changes the amino acid at the specified position
    # to a different random amino acid

    seq_str = str(sequence)

    if position < 0 or position >= len(seq_str):
        raise ValueError("Position is out of range for this sequence.")

    current_aa = seq_str[position]

    # Ensure it mutates to a different amino acid
    possible_mutations = [aa for aa in AMINO_ACIDS if aa != current_aa]
    new_aa = random.choice(possible_mutations)

    # Replace the amino acid at the specified position
    mutated_list = list(seq_str)
    mutated_list[position] = new_aa

    return Seq("".join(mutated_list))


def random_mutation(sequence):
    # Changes a random amino acid to a different amino acid

    seq_str = str(sequence)

    if len(seq_str) == 0:
        raise ValueError("Sequence is empty.")

    random_pos = random.randint(0, len(seq_str) - 1)

    print(f"Mutating position {random_pos} (0-indexed)...")

    return mutate_position(seq_str, random_pos)


def run_mutations(sequence, iterations=10000):

    # Create a list of zeros, one for each position
    mutation_counts = [0] * len(sequence)

    for i in range(iterations):

        # Start each iteration with the original sequence
        current_seq = sequence

        # Go through each amino acid position sequentially
        for position in range(len(sequence)):

            # Roll a 100-sided die
            dice_roll = random.randint(1, 100)

            # If the die lands on 1, mutate
            if dice_roll == 1:

                current_seq = mutate_position(current_seq, position)

                # Count how many times this position mutated
                mutation_counts[position] += 1

    return mutation_counts


if __name__ == "__main__":

    original_seq = Seq(
        "PAGFRQVCKAFMAEPMFNTATQHEYHCCTPTRMYINTPVNRAFGCFDGPIIFGWPPEREIYPNQTS"
        "FAHGFWVHQIWLGFNGNVIFYMCADPSKCIYMYVARQMGPDKLYVGMEWVGPPESMSMAHAWVTHC"
        "VWCLCSSDLMGVKFMRFHYPKMDFHNFFCPFWVCFQFDWNLSTECKKLRQYCKHVFEKKGFWKVNC"
        "YTYWTELFKCIYGSLKYYNREGYNMRCNMVEYCENQFDVFNPKNMWFRFIQGHSCLYYCFGWVVIE"
        "TRPYQATGQWDDDNGPRFPGLWTCKGASPGRNKACVWMTNDDHVSVVFFDSVKPEWVYTVARHWTV"
        "TGSDVVFTLHSFAEHVYGPIYMWDYEYGDKRHPMINGRAWMRRQICQHHTRHHKKGDKWQWPHMAG"
        "KNQLKFLWSECMPHTSQFEAFLTCKPEQLFQASSDKEYMCFTLTQYNYVPRMTMVRLGFPSMPGCG"
        "GHMDFQRCKMKFKWNQVSIAISVYKRMPFFKAHQKELDPAINAPDLCIFHSISTGCAMYHQSQSWE"
    )

    iterations = 100000

    # Run simulation
    results = run_mutations(original_seq, iterations)

    # Calculate statistics
    mean = statistics.mean(results)
    std_dev = statistics.stdev(results)

    print(f"Mean mutations per position: {mean:.2f}")
    print(f"Standard deviation: {std_dev:.2f}")

    # Print mutation counts for each amino acid position
    print("\nPosition | Amino Acid | Mutation Count")
    print("--------------------------------------")

    for pos, count in enumerate(results):
        print(f"{pos + 1:<8} | {original_seq[pos]:<10} | {count}")

    print(f"\nTotal mutations: {sum(results)}")

    # Plot mutation frequencies
    positions = list(range(1, len(original_seq) + 1))

    plt.figure(figsize=(12, 5))

    plt.plot(positions, results, marker="o", markersize=3)

    # Expected number of mutations per position
    expected = iterations / 100

    plt.axhline(
        y=expected,
        color="red",
        linestyle="--",
        label=f"Expected mutations ({expected:.0f})",
    )

    plt.xlabel("Amino Acid Position")
    plt.ylabel("Times Mutated")
    plt.title(f"Mutation Count at Each Amino Acid Position ({iterations:,} iterations)")

    plt.xticks(range(1, len(original_seq) + 1, 5))
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.show()
