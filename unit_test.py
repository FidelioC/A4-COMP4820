import de_bruijn


def test_possible_kmers():
    print(de_bruijn.create_possible_kmers(5, "AGCCGATCAT"))


def test_one_graph_read():
    possible_kmers = de_bruijn.create_possible_kmers(5, "AGCCGATCAT")
    graph = de_bruijn.create_graph(possible_kmers)
    print(graph)
    de_bruijn.visualize_graph(graph)


def test_all_graph_read():
    all_reads = [
        "AGCCGATCAT",
        "CCGATCATC",
        "GCCGAT",
        "AGCCGA",
        "GATCATGGTCA",
        "GATCAAGGTCA",
    ]
    graph = de_bruijn.create_graph_all_reads(all_reads, 5)
    de_bruijn.visualize_graph(graph)


def main():
    # test_possible_kmers()
    # test_one_graph_read()
    test_all_graph_read()


if __name__ == "__main__":
    main()
