import de_bruijn


def test_possible_kmers():
    print(de_bruijn.create_possible_kmers(5, "CCGATCATC"))


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
    graph = de_bruijn.create_graph_list_reads(all_reads, 5)
    # de_bruijn.visualize_graph(graph)
    print(graph)
    out_degree = de_bruijn.find_all_outdegree_zero(graph)
    in_degree = de_bruijn.find_all_indegree_zero(graph)
    print(f"total nodes: {len(graph)}")
    print(f"out degree: {len(out_degree)}")
    print(f"in degree: {len(in_degree)}")

    de_bruijn.visualize_graph(graph, "./graph_output/class_example_output")


def test_degree_zero():
    all_reads = [
        "AGCCGATCAT",
        "CCGATCATC",
        "GCCGAT",
        "AGCCGA",
        "GATCATGGTCA",
        "GATCAAGGTCA",
    ]
    graph = de_bruijn.create_graph_list_reads(all_reads, 5)
    out_degree = de_bruijn.find_all_outdegree_zero(graph)
    in_degree = de_bruijn.find_all_indegree_zero(graph)
    print(f"graph: {graph}")
    print(f"out degree: {out_degree}")
    print(f"in degree: {in_degree}")


def test_graph_progress():
    file_name = "./simulated-reads/covid.fq"
    k = 37
    graph = de_bruijn.create_graph_file_reads(file_name, k)

    out_degree = de_bruijn.find_all_outdegree_zero(graph)
    in_degree = de_bruijn.find_all_indegree_zero(graph)
    tip_outdegree, tip_indegree = de_bruijn.find_tip(graph, k)

    print(f"total nodes: {len(graph)}")
    print(f"out degree: {len(out_degree)}")
    print(f"in degree: {len(in_degree)}")
    # print(f"graph: {graph}")
    # print(f"tip outdegree, tip indegree {tip_outdegree,tip_indegree}")
    # graph, total_prune_tips, total_bubbles_removed = de_bruijn.remove_all_errors(
    #     graph, k
    # )
    # print(f"pruned {total_prune_tips} nodes")
    # print(f"removed {total_bubbles_removed} nodes in bubbles")
    # de_bruijn.visualize_graph(graph, "./graph_output/graph_output_covidfq")


def test_all_graph_errors():
    all_reads = [
        "AGCCGATCAT",
        "CCGATCATC",
        "GCCGAT",
        "AGCCGA",
        "GATCATGGTCA",
        "GATCAAGGTCA",
    ]
    graph = de_bruijn.create_graph_list_reads(all_reads, 5)
    # de_bruijn.visualize_graph(graph)
    # print(graph)
    out_degree = de_bruijn.find_all_outdegree_zero(graph)
    in_degree = de_bruijn.find_all_indegree_zero(graph)
    print(f"total nodes: {len(graph)}")
    print(f"out degree: {len(out_degree)}")
    print(f"in degree: {len(in_degree)}")

    de_bruijn.visualize_graph(graph, "./graph_output/class_example_output")

    de_bruijn.find_tip(graph)


def test_possible_inoutgoing():
    possible_outgoing = de_bruijn.generate_possible_outgoing("ATCAT")
    print(f"possible_outgoing: {possible_outgoing}")
    possible_ingoing = de_bruijn.generate_possible_ingoing("ATCAT")
    print(f"possible_ingoing:{possible_ingoing}")


def test_find_parent():
    all_reads = [
        "AGCCGATCAT",
        "CCGATCATC",
        "GCCGAT",
        "AGCCGA",
        "GATCATGGTCA",
        "GATCAAGGTCA",
    ]
    graph = de_bruijn.create_graph_list_reads(all_reads, 5)

    # print(de_bruijn.find_parent(graph, "GGTCA"))
    all_nodes_indegree = de_bruijn.find_nodes_indegree(graph)
    # print(all_nodes_indegree)
    print(f"ATCAT: {de_bruijn.find_tip_traverse_parents(graph, 'ATCAT', 1, 5, [])}")
    print(f"TCATC: {de_bruijn.find_tip_traverse_parents(graph, 'TCATC', 1, 5, [])}")
    print(
        de_bruijn.find_tip_traverse_children(
            graph, "AGCCG", 1, 5, [], all_nodes_indegree
        )
    )
    de_bruijn.visualize_graph(graph, "./test_graph_tip")


def test_find_tips():
    all_reads = ["AAGCCGATCAT", "CTCGGATCA", "ATCATCG", "ATCATGGTCA"]

    graph = de_bruijn.create_graph_list_reads(all_reads, 5)
    print(de_bruijn.find_tip(graph, 5))

    de_bruijn.visualize_graph(graph, "./test_graph_tip")


def test_remove_all_tips():
    all_read_tips = ["AAGCCGATCAT", "CTCGGATCA", "ATCATCG", "ATCATGGTCA"]
    graph = de_bruijn.create_graph_list_reads(all_read_tips, 5)
    graph, total_prune_outdegree, total_prune_indegree = (
        de_bruijn.remove_graph_all_tips(graph, 5)
    )

    print(f"pruned {total_prune_outdegree + total_prune_indegree} nodes")

    de_bruijn.visualize_graph(graph, "./test_graph_tip")


def test_out_degree_greater_zero():
    all_read_tips = [
        "AAGCCGATCAT",
        "CTCGGATCA",
        "ATCATCGGTCA",
        "ATCATGGTCA",
        "GATCAAGGTCA",
    ]
    graph = de_bruijn.create_graph_list_reads(all_read_tips, 5)

    print(de_bruijn.find_all_outdegree_greater_zero(graph))


def test_find_bubbles():
    all_read_tips = [
        "AAGCCGATCAT",
        "CTCGGATCA",
        "ATCATCGGTCA",
        "ATCATGGTCA",
        "GATCAAGGTCA",
        "GGTCACG",
    ]
    graph = de_bruijn.create_graph_list_reads(all_read_tips, 5)
    # print(graph)
    # graph, total_prune_outdegree, total_prune_indegree = (
    #     de_bruijn.remove_graph_all_tips(graph, 5)
    # )

    # print(f"pruned {total_prune_outdegree + total_prune_indegree} nodes")

    print(de_bruijn.traverse_find_bubble_iterative(graph, "ATCAT"))  # has bubbles
    print(de_bruijn.traverse_find_bubble_iterative(graph, "GATCA"))  # has bubbles
    print(de_bruijn.traverse_find_bubble_iterative(graph, "AAGCC"))  # no bubble
    print("\n")
    print(f"all bubbles: {de_bruijn.find_bubble(graph)}")
    de_bruijn.visualize_graph(graph, "./test_graph_tip")


def test_find_lowest_weight_children():
    children = {"ATCAA": 1, "ATCAT": 3, "TTTTT": 5}
    max_weight = de_bruijn.find_lowest_weight_children(children)
    print(max_weight)


def test_remove_node_graph():
    all_read_tips = [
        "AAGCCGATCAT",
        "CTCGGATCA",
        "ATCATCGGTCA",
        "ATCATGGTCA",
        "GATCAAGGTCA",  # bubble GATCA
        "GGTCACG",
        "GATCAT",  # this will make GATCA -> ATCAT weight = 2
        "GATCATG",  # this will make GATCA -> ATCAT weight = 3 and ATCAT -> TCATG weight = 2
        "GATCAGG",  # this will make GATCA -> ATCAG weight = 1
        "GATCAA",  # will make GATCA -> ATCAA weight = 2
    ]
    graph = de_bruijn.create_graph_list_reads(all_read_tips, 5)

    graph = de_bruijn.remove_node_from_graph("GATCA", graph)
    de_bruijn.visualize_graph(graph, "./test_graph")


def test_remove_bubbles():
    all_read_tips = [
        "AAGCCGATCAT",
        "CTCGGATCA",
        "ATCATCGGTCA",
        "ATCATGGTCA",
        "GATCAAGGTCA",  # bubble GATCA
        "GGTCACG",
        "GATCAT",  # this will make GATCA -> ATCAT weight = 2
        "GATCATG",  # this will make GATCA -> ATCAT weight = 3 and ATCAT -> TCATG weight = 2
        "GATCAGG",  # this will make GATCA -> ATCAG weight = 1
        "GATCAA",  # will make GATCA -> ATCAA weight = 2
    ]
    graph = de_bruijn.create_graph_list_reads(all_read_tips, 5)

    print(de_bruijn.remove_bubble(graph, 5))

    de_bruijn.visualize_graph(graph, "./test_graph")


def test_remove_all_errors():
    all_read_tips = [
        "GTGTAACATTAGGGAGGACTTGAAAGAGCCACCACATATTTCACCGAGGCCACGCGGAGTACGATCGAGTGTACA",
        "GTGTAACATTAGGGAGGACTTGAAAGAGCCACCACATTTTCACCGAGGCCACGCGGAGTACGATCGAGTGTACA",
        "GTGTAACATTAGGGAGGACTTGAAAGAGCCACCACATTTTCACCGAGGCCACGCGGAGTACGATCGAGTGTACA",
    ]
    k_file = 37
    # k_list = 5
    file_name = "./simulated-reads/covid.fasta_simulated-errors-tips.fq"
    # graph = de_bruijn.create_graph_list_reads(all_read_tips, k_file)
    graph = de_bruijn.create_graph_file_reads(file_name, k_file)
    # de_bruijn.visualize_graph(graph, "./graph_before_error_removal")
    out_degree = de_bruijn.find_all_outdegree_zero(graph)
    in_degree = de_bruijn.find_all_indegree_zero(graph)

    print(f"total nodes: {len(graph)}")
    print(f"out degree: {len(out_degree)}")
    print(f"in degree: {len(in_degree)}")

    # graph, indegree, outdegree = de_bruijn.remove_graph_all_tips(graph, k_file)
    # print(indegree + outdegree)

    graph, total_tips_pruned, total_bubbles_removed = de_bruijn.remove_all_errors(
        graph, k_file
    )

    print(
        f"final all_outdegree greater zero: {de_bruijn.find_all_outdegree_greater_zero(graph)}"
    )
    print(
        f"total_tips_pruned {total_tips_pruned}, total_bubbles_removed {total_bubbles_removed}"
    )
    de_bruijn.visualize_graph(graph, "./test_graph")


def main():
    # test_possible_kmers()
    # test_one_graph_read()
    # test_all_graph_read()
    # test_degree_zero()
    # test_graph_progress()
    # test_all_graph_errors()
    # test_possible_inoutgoing()
    # test_find_tips()
    # test_find_parent()
    # test_remove_all_tips()
    # test_out_degree_greater_zero()
    # test_find_bubbles()
    # test_find_lowest_weight_children()
    # test_remove_node_graph()
    # test_remove_bubbles()
    test_remove_all_errors()


if __name__ == "__main__":
    main()
