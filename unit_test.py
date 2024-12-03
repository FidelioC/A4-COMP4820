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
    graph = de_bruijn.create_graph_all_reads(all_reads, 5)
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
    graph = de_bruijn.create_graph_all_reads(all_reads, 5)
    out_degree = de_bruijn.find_all_outdegree_zero(graph)
    in_degree = de_bruijn.find_all_indegree_zero(graph)
    print(f"graph: {graph}")
    print(f"out degree: {out_degree}")
    print(f"in degree: {in_degree}")


def test_graph_progress():
    file_name = "./simulated-reads/covid.fasta_simulated-errors-tips.fq"
    k = 37
    graph = de_bruijn.progress(file_name, k)

    out_degree = de_bruijn.find_all_outdegree_zero(graph)
    in_degree = de_bruijn.find_all_indegree_zero(graph)
    tip_outdegree, tip_indegree = de_bruijn.find_tip(graph, k)

    print(f"total nodes: {len(graph)}")
    print(f"out degree: {len(out_degree)}")
    print(f"in degree: {len(in_degree)}")
    # print(f"graph: {graph}")
    # print(f"tip outdegree, tip indegree {tip_outdegree,tip_indegree}")
    graph, total_prune_outdegree, total_prune_indegree = (
        de_bruijn.remove_graph_all_tips(graph, k)
    )
    print(f"pruned {total_prune_outdegree + total_prune_indegree} nodes")
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
    graph = de_bruijn.create_graph_all_reads(all_reads, 5)
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


def test_find_tips():
    all_reads = [
        "AGCCGATCAT",
        "CCGATCATC",
        "GCCGAT",
        "AGCCGA",
        "GATCATGGTCA",
        "GATCAAGGTCA",
    ]
    graph = de_bruijn.create_graph_all_reads(all_reads, 5)
    print(de_bruijn.find_tip(graph, 5))

    de_bruijn.visualize_graph(graph, "./test_graph_tip")


def test_find_parent():
    all_reads = [
        "AGCCGATCAT",
        "CCGATCATC",
        "GCCGAT",
        "AGCCGA",
        "GATCATGGTCA",
        "GATCAAGGTCA",
    ]
    graph = de_bruijn.create_graph_all_reads(all_reads, 5)

    # print(de_bruijn.find_parent(graph, "GGTCA"))
    all_nodes_indegree = de_bruijn.find_nodes_indegree(graph)
    print(all_nodes_indegree)
    print(de_bruijn.find_tip_traverse_parents(graph, "ATCAT", 1, 5, []))
    print(
        de_bruijn.find_tip_traverse_children(
            graph, "AGCCG", 1, 5, [], all_nodes_indegree
        )
    )


def test_remove_all_tips():
    all_read_tips = ["AAGCCGATCAT", "CTCGGATCA", "ATCATCG", "ATCATGGTCA"]
    graph = de_bruijn.create_graph_all_reads(all_read_tips, 5)
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
    graph = de_bruijn.create_graph_all_reads(all_read_tips, 5)

    print(de_bruijn.find_all_outdegree_greater_zero(graph))


def test_remove_bubbles():
    all_read_tips = [
        "AAGCCGATCAT",
        "CTCGGATCA",
        "ATCATCGGTCA",
        "ATCATGGTCA",
        "GATCAAGGTCA",
    ]
    graph = de_bruijn.create_graph_all_reads(all_read_tips, 5)
    # print(graph)
    # graph, total_prune_outdegree, total_prune_indegree = (
    #     de_bruijn.remove_graph_all_tips(graph, 5)
    # )

    # print(f"pruned {total_prune_outdegree + total_prune_indegree} nodes")

    print(de_bruijn.traverse_find_bubble(graph, "GATCA"))
    # print(de_bruijn.find_bubble(graph))
    de_bruijn.visualize_graph(graph, "./test_graph_tip")


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
    test_remove_bubbles()


if __name__ == "__main__":
    main()
