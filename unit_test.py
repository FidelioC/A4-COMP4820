import de_bruijn


def main():
    # print(de_bruijn.create_possible_kmers(5, "AGCCGATCAT"))
    possible_kmers = de_bruijn.create_possible_kmers(4, "ATCGATCGA")
    graph = de_bruijn.create_graph(possible_kmers)
    print(graph)

    de_bruijn.visualize_graph(graph)


if __name__ == "__main__":
    main()
