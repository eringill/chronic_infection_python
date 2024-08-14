import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import functions


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="""SMDP: SARS-CoV-2 Mutation Distribution Profiler for rapid estimation of mutational histories.
Developed by the CAMEO team of CoVaRR-Net.

For the web app, visit: https://eringill.shinyapps.io/covid_mutation_distributions/
Source code: https://github.com/eringill/chronic_infection_python
Please cite: https://doi.org/10.48550/arXiv.2407.11201""",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "mutations",
        help="Comma-separated list of mutations or path to a file containing mutations",
    )
    parser.add_argument(
        "--bin-size",
        choices=["genes_split", "gene", "500", "1000"],
        default="gene",
        help="Bin size for analysis (default: gene)",
    )
    parser.add_argument(
        "--output",
        choices=["text", "json"],
        default="text",
        help="Output format (default: text)",
    )
    parser.add_argument(
        "--plot", action="store_true", help="Generate a plot of mutation distribution"
    )
    parser.add_argument(
        "--plot-output",
        default="mutation_distribution.png",
        help="Output file for the plot (default: mutation_distribution.png)",
    )
    parser.add_argument(
        "--color-palette",
        choices=["plasma", "viridis", "inferno", "seaborn"],
        default="plasma",
        help="Color palette for the plot (default: plasma)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed information during analysis",
    )
    return parser.parse_args()


def load_mutations(mutations_input: str) -> List[str]:
    is_file = False
    try:
        is_file = Path(mutations_input).is_file()
    except OSError as e:
        print(f"Error: Unable to read mutations as file: {e}")
    if is_file:
        with open(mutations_input, "r") as f:
            mutations = f.read().strip()
    else:
        mutations = mutations_input
    return functions.parse_user_input(mutations)


def load_distribution_data() -> Dict[str, Tuple[List[int], int]]:
    data_dir = Path(__file__).parent / "data"
    distributions = {
        "chronic": "chronicnucl.tsv",
        "deer": "deernucl.tsv",
        "global": "globalnucl.tsv",
        "global_late": "globallatenucl.tsv",
    }

    data = {}
    for name, filename in distributions.items():
        file_path = data_dir / filename
        data[name], data[f"total_{name}"] = functions.parse_mutation_files(file_path)

    return data


def calculate_likelihoods(
    bin_size: str,
    mutations: List[str],
    distribution_data: Dict[str, Tuple[List[int], int]],
) -> Tuple[List[Tuple[float, str]], str]:
    likelihood_list, most_likely = functions.most_likely(
        bin_size,
        distribution_data["global"],
        distribution_data["global_late"],
        distribution_data["chronic"],
        distribution_data["deer"],
        mutations,
    )
    return likelihood_list, most_likely


def analyze_mutations(
    mutations: str, bin_size: str, verbose: bool
) -> Tuple[Dict[str, any], List[str], Dict[str, Tuple[List[int], int]]]:
    mut_list = load_mutations(mutations)
    transitions, transversions = functions.transition_or_transversion(mutations)

    if verbose:
        print(f"Analyzing {len(mut_list)} mutations...")
        print(f"Transitions: {transitions}, Transversions: {transversions}")

    distribution_data = load_distribution_data()
    likelihood_list, most_likely = calculate_likelihoods(
        bin_size, mutations, distribution_data
    )

    results = {
        "mutations_count": len(mut_list),
        "transition_transversion_ratio": transitions / transversions
        if transversions
        else None,
        "likelihoods": {
            dist.replace("_", " "): likelihood[0]
            for dist, likelihood in zip(distribution_data.keys(), likelihood_list)
        },
        "best_fit": most_likely,
        "times_more_likely": functions.times_more_likely(likelihood_list)[0],
        "compared_to": functions.times_more_likely(likelihood_list)[1].replace(
            "_", " "
        ),
        "mutator_lineage": functions.mut_lineage_parsing(mutations),
    }

    if verbose:
        print("Analysis complete.")

    return results, mut_list, distribution_data


def print_results(results: Dict[str, any]) -> None:
    print(f"Number of mutations: {results['mutations_count']}")
    print(
        f"Transition/Transversion ratio: {results['transition_transversion_ratio']:.2f}"
    )
    print("\nLog Likelihoods:")
    for dist, likelihood in results["likelihoods"].items():
        print(f"  {dist}: {likelihood:.2f}")
    print(f"\nBest fit distribution: {results['best_fit']}")
    print(
        f"({results['times_more_likely']:.2f} times more likely than the {results['compared_to']} distribution)"
    )

    print("\nMutator lineage analysis:")
    if results["mutator_lineage"][0]:
        print(f"  Confirmed: {results['mutator_lineage'][0]}")
    if results["mutator_lineage"][1]:
        print(f"  Potential: {results['mutator_lineage'][1]}")
    if not any(results["mutator_lineage"]):
        print("  No mutator lineage detected")


def generate_plot(
    mut_list: List[str],
    bin_size: str,
    distribution_data: Dict[str, Tuple[List[int], int]],
    color_palette: str,
    output_file: str,
) -> None:
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print(
            "Error: matplotlib and seaborn are required for plotting. Please install them and try again."
        )
        return

    counts, bins = functions.make_bins(mut_list, bin_size)

    plt.figure(figsize=(12, 6))
    sns.set_style("whitegrid")
    sns.set_palette(color_palette)

    for name, data in distribution_data.items():
        if not name.startswith("total_"):
            dist_counts, _ = functions.make_bins(data, bin_size)
            total = distribution_data[f"total_{name}"]
            plt.plot(
                bins, [x / total for x in dist_counts], label=name.replace("_", " ")
            )

    plt.plot(
        bins, [x / len(mut_list) for x in counts], label="input mutations", linewidth=2
    )

    plt.xlabel("Genome Position")
    plt.ylabel("Proportion of Mutations")
    plt.title("Distribution of Mutations Across Genome")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Plot saved as {output_file}")


def main() -> None:
    args = parse_arguments()

    if args.plot and args.bin_size not in ["500", "1000"]:
        print(
            "Error: Plotting is only available for integer bin sizes. Set `--bin-size` to 500 or 1000 to plot."
        )
        exit(1)

    results, mut_list, distribution_data = analyze_mutations(
        args.mutations, args.bin_size, args.verbose
    )

    if args.output == "text":
        print_results(results)
    elif args.output == "json":
        print(json.dumps(results, indent=2))

    if args.plot:
        generate_plot(
            mut_list,
            args.bin_size,
            distribution_data,
            args.color_palette,
            args.plot_output,
        )


if __name__ == "__main__":
    main()
