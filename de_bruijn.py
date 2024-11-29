import graphviz

"""'
graph = {
    AGCCG: 
    GCCGA:


}
"""


def find_all_outdegree_zero():
    """
    find all nodes that have outdegree of zero
    """


def create_possible_kmers(k, read):
    """generate all possible k-mers from one read"""
    all_kmers = []
    for i in range(len(read) - k + 1):
        all_kmers.append(read[i : i + k])

    return all_kmers


def check_node_in_graph(graph: dict, kmer_1, kmer_2):
    """check if Vi or Vj is in the graph"""
    if kmer_1 in graph and kmer_2 in graph:
        return (kmer_1, kmer_2)  # both exist in graph
    elif kmer_1 in graph:
        return (kmer_1, None)  # kmer_1 exist in graph but kmer_2 doesn't exist
    elif kmer_2 in graph:
        return (None, kmer_2)  # kmer_2 exist in graph but kmer_1 doesn't exist
    else:
        return (None, None)  # both doesn't exist in graph


def edit_graph_both_exist(graph, kmer_1, kmer_2):
    if kmer_2 not in graph[kmer_1]:
        graph[kmer_1][kmer_2] = 1
    else:
        graph[kmer_1][kmer_2] += 1


def edit_graph_kmer1_exist(graph, kmer_1, kmer_2):
    if kmer_2 not in graph[kmer_1]:
        graph[kmer_1][kmer_2] = 1
    else:
        graph[kmer_1][kmer_2] += 1
    graph[kmer_2] = {}


def edit_graph_kmer2_exist(graph, kmer_1, kmer_2):
    if kmer_1 not in graph[kmer_2]:
        graph[kmer_2][kmer_1] = 1
    else:
        graph[kmer_2][kmer_1] += 1
    graph[kmer_1] = {}


def edit_graph_both_dont_exist(graph, kmer_1, kmer_2):
    graph[kmer_1] = {}
    graph[kmer_1][kmer_2] = 1
    graph[kmer_2] = {}


def create_graph(curr_graph, all_kmers):
    """
    STEPS:
    1) Take contiguos pairs of k-mers. Vi, Vj, Eij
    2) check if Vi or Vj is in the graph
        - if both: add 1 to Eij.
        - if one of them: add either Vi or Vj to G.
        - if none: add Vi, Vj, Eij to G.
    """
    graph = curr_graph
    index = 0
    check_result = None
    while index + 1 < len(all_kmers):
        kmer_1 = all_kmers[index]
        kmer_2 = all_kmers[index + 1]
        check_result = check_node_in_graph(graph, kmer_1, kmer_2)
        print(check_result)
        if check_result == (kmer_1, kmer_2):
            # print("Both exist")
            edit_graph_both_exist(graph, kmer_1, kmer_2)
        elif check_result == (kmer_1, None):
            # print("kmer1 exist")
            edit_graph_kmer1_exist(graph, kmer_1, kmer_2)
        elif check_result == (None, kmer_2):
            # print("kmer2 exist")
            edit_graph_kmer2_exist(graph, kmer_1, kmer_2)
        else:
            # print("none exist")
            edit_graph_both_dont_exist(graph, kmer_1, kmer_2)

        print(f"kmer1: {kmer_1}, kmer2: {kmer_2}, graph: {graph}\n")
        index += 1

    return graph


def create_graph_all_reads(all_reads, k):
    curr_graph = {}
    for read in all_reads:
        print(f"============= CURRENT READ {read} ============= \n")
        print(f"graph: {curr_graph}")
        kmers = create_possible_kmers(k, read)
        curr_graph = create_graph(curr_graph, kmers)
        print(f"KMERS: {kmers}")
    return curr_graph


def visualize_graph(graph: dict):
    dot = graphviz.Digraph(format="svg")

    for node, neighbours in graph.items():
        dot.node(node)
        for neighbour in neighbours:
            weight = neighbours[neighbour]
            dot.node(neighbour)
            dot.edge(node, neighbour, label=str(weight))

    dot.render("de_bruijn_graph_output")


def main():
    print(create_possible_kmers(5, "AGCCGATCAT"))


if __name__ == "__main__":
    main()
