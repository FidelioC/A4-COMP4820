import graphviz
import sys, time

from rich.progress import Progress
from Bio import SeqIO

ALL_POSSIBLE_CHARS = ["A", "C", "G", "T"]


# ====================== graph construction code ================== "
def find_all_outdegree_zero(graph):
    """
    find all nodes that have outdegree of zero
    """
    outdegree_zero = []
    for node, neighbours in graph.items():
        if neighbours == {}:
            outdegree_zero.append(node)

    return outdegree_zero


def find_all_indegree_zero(graph):
    """
    find all nodes that have indegree of zero
    """
    init_nodes = find_nodes_indegree(graph)
    indegree_zero = []

    for node, value in init_nodes.items():
        if value == 0:
            indegree_zero.append(node)

    # print(f"init_nodes {init_nodes}")

    return indegree_zero


def find_nodes_indegree(graph):
    init_nodes = {node: 0 for node in graph}  # start all nodes with indegree 0

    for node, neighbours in graph.items():
        for destination in neighbours:
            if destination in init_nodes:
                init_nodes[destination] += 1

    return init_nodes


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
        # print(check_result)
        if check_result == (kmer_1, kmer_2):
            # print("Both exist")
            edit_graph_both_exist(graph, kmer_1, kmer_2)
            # weight = graph[kmer_1][kmer_2]
        elif check_result == (kmer_1, None):
            # print("kmer1 exist")
            edit_graph_kmer1_exist(graph, kmer_1, kmer_2)
            # dot.node(kmer_2)
            # weight = 1
        elif check_result == (None, kmer_2):
            # print("kmer2 exist")
            edit_graph_kmer2_exist(graph, kmer_1, kmer_2)
            # weight = 1
        else:
            # print("none exist")
            edit_graph_both_dont_exist(graph, kmer_1, kmer_2)
            # dot.node(kmer_1)
            # dot.node(kmer_2)
            # weight = 1

        # print(f"kmer1: {kmer_1}, kmer2: {kmer_2}, graph: {graph}\n")
        index += 1

        # dot.edge(kmer_1, kmer_2, label=str(weight))
        # print(graph)

    return graph


def create_graph_all_reads(all_reads, k):
    curr_graph = {}
    # dot = graphviz.Digraph(format="svg")
    for read in all_reads:
        # print(f"============= CURRENT READ {read} ============= \n")
        # print(f"graph: {curr_graph}")
        kmers = create_possible_kmers(k, read)
        curr_graph = create_graph(curr_graph, kmers)
        # print(f"KMERS: {kmers}")
    # dot.render("de_bruijn_graph_output")
    return curr_graph


def visualize_graph(graph: dict, output_file):
    dot = graphviz.Digraph(format="svg")

    for node, neighbours in graph.items():
        dot.node(node)
        for neighbour in neighbours:
            weight = neighbours[neighbour]
            dot.node(neighbour)
            dot.edge(node, neighbour, label=str(weight))

    # print("finished creating visualized graph, waiting for render...")
    dot.render(output_file)


# =============== errors tip ====================
def find_tip(graph, k):
    """find all tip nodes in the graph"""
    nodes_indegree_zero = find_all_indegree_zero(graph)
    nodes_outdegree_zero = find_all_outdegree_zero(graph)

    all_potential_tips = nodes_indegree_zero + nodes_outdegree_zero
    # print(f"nodes_indegree_zero {nodes_indegree_zero}")
    # print(f"nodes_outdegree_zero {nodes_outdegree_zero}")
    # print(f"all potential tips: {all_potential_tips}")

    tip_outdegree = find_tip_outdegree_zero(graph, nodes_outdegree_zero, k)
    tip_indegree = find_tip_indegree_zero(graph, nodes_indegree_zero, k)
    # print(f"tip outdegree zero: {tip_outdegree}")
    # print(f"tip indegree zero: {tip_indegree}")
    return tip_outdegree, tip_indegree


def remove_graph_all_tips(graph, k):
    """remove tips from the graph"""
    tip_outdegree, tip_indegree = find_tip(graph, k)

    total_prune_outdegree = len(tip_outdegree)
    total_prune_indegree = len(tip_indegree)
    total_tips = total_prune_outdegree + total_prune_indegree

    while total_tips > 0:  # while tip still exist in the graph
        # print(f"tips: {tip_outdegree, tip_indegree}")
        graph = remove_tip_indegree(tip_indegree, graph)
        graph = remove_tip_outdegree(tip_outdegree, graph)
        # visualize_graph(graph, "./test_graph_tip")
        tip_outdegree, tip_indegree = find_tip(
            graph, k
        )  # find another possible tips in the graph

        total_prune_outdegree += len(tip_outdegree)
        total_prune_indegree += len(tip_indegree)
        total_tips = len(tip_outdegree) + len(tip_indegree)

    return graph, total_prune_outdegree, total_prune_indegree


def remove_tip_outdegree(tip_outdegree, graph):
    """remove tips that has outdegree 0
    (we will have to edit the parent's node, where the edge is pointing to this node)
    """
    for tip in tip_outdegree:
        # print(f"removing tip outdegree {tip}")
        parents = find_parent(graph, tip)
        for parent in parents:
            del graph[parent][tip]  # delete the parent's edge pointing to the tip
        del graph[tip]  # delete the tip itself

    return graph


def remove_tip_indegree(tip_indegree, graph):
    """remove tips that has indegree 0"""
    for tip in tip_indegree:
        # print(f"removing tip indegree {tip}")
        del graph[tip]

    return graph


def find_tip_indegree_zero(graph, nodes_indegree_zero, k):
    """find indegree tips"""
    tip_path = []  # might be used later
    tip_nodes = []
    all_nodes_indegree = find_nodes_indegree(graph)
    for node in nodes_indegree_zero:
        tip_path = []
        tip_path = find_tip_traverse_children(
            graph, node, 1, k, tip_path, all_nodes_indegree
        )
        if tip_path:
            tip_nodes.append(node)
    return tip_nodes


def find_tip_traverse_children(graph, node, depth, k, path: list, all_nodes_indegree):
    """check if node with in degree zero is a tip"""
    path.append(node)

    if depth >= k:
        return False

    children = list(graph.get(node, {}).keys())  # Get the children of the current node
    for child in children:
        if all_nodes_indegree[child] > 1:  # If a child has indegree > 1
            path.append(child)
            return path
        # Recursively check the child's children
        result = find_tip_traverse_children(
            graph, child, depth + 1, k, path, all_nodes_indegree
        )
        if result:
            return result

    return False


def find_tip_outdegree_zero(graph, nodes_outdegree_zero, k):
    """find out degree tips"""
    tip_path = []  # might be used later
    tip_nodes = []
    for node in nodes_outdegree_zero:
        tip_path = []
        tip_path = find_tip_traverse_parents(graph, node, 1, k, tip_path)
        if tip_path:
            tip_nodes.append(node)
    return tip_nodes


def find_tip_traverse_parents(graph, node, depth, k, path: list):
    """check if node with out degree zero is a tip"""
    path.append(node)
    if depth >= k:
        return False

    parents = find_parent(graph, node)
    # print(f"parents: {parents}")
    for parent in parents:
        curr_parent_outdegree = len(graph[parent])
        if curr_parent_outdegree > 1:
            path.append(parent)
            return path

        result = find_tip_traverse_parents(graph, parent, depth + 1, k, path)
        if result:
            return result

    return False


def generate_possible_outgoing(node):
    """generate all possible outgoing node given a node"""
    prune_first_char = node[1:]
    all_possible_outgoing = []
    for char in ALL_POSSIBLE_CHARS:
        all_possible_outgoing.append(prune_first_char + char)

    return all_possible_outgoing


def generate_possible_ingoing(node):
    """generate all possible ingoing node given a node"""
    prune_last_char = node[:-1]
    all_possible_ingoing = []
    for char in ALL_POSSIBLE_CHARS:
        all_possible_ingoing.append(char + prune_last_char)

    return all_possible_ingoing


def find_parent(graph, target_node):
    all_possible_ingoing = generate_possible_ingoing(target_node)
    all_parents = []
    for node in all_possible_ingoing:
        if node in graph and target_node in graph[node]:
            all_parents.append(node)

    return all_parents


# ============= errors bubbles =============
def remove_bubble(graph):
    potential_bubbles = find_bubble(graph)


def find_bubble(graph):
    """find potential bubbles in nodes with outdegree > 1"""
    potential_bubbles = {}
    outdegree_greater_zero = find_all_outdegree_greater_zero(graph)
    for node in outdegree_greater_zero:
        result = traverse_find_bubble(graph, node)
        potential_bubbles[node] = result

    return potential_bubbles


def traverse_find_bubble(graph, curr_node, path=[], paths_end={}):
    """given a node (outdegree > 1), will traverse the graph,
    see if the path ends in the same node"""
    path.append(curr_node)

    curr_node_children = graph[curr_node]

    for child in curr_node_children:
        if child in paths_end:  # if have reached common nodes
            # append the previous path with the current child
            paths_end[child].append(path + [child])
            if len(paths_end[child]) > 1:
                return paths_end[child]

        else:  # else, haven't reached common nodes, just keep recording the path for now
            paths_end[child] = [path + [child]]
        # print(paths_end)
        result = traverse_find_bubble(graph, child, path[:], paths_end)
        if result:
            return result

    return None


def find_all_outdegree_greater_zero(graph):
    """
    find all nodes that have outdegree greater than zero
    """
    outdegree_greater_zero = []
    for node, neighbours in graph.items():
        if len(neighbours) > 1:
            outdegree_greater_zero.append(node)

    return outdegree_greater_zero


# =============== contigs ===================


# =============== progress ==================
def progress(read_file, k):
    curr_graph = {}
    with open(read_file) as reads:
        # Count all records in the file so we know to track
        # overall progress (FASTQ has 4 lines per record).
        num_reads = sum(1 for _ in open(read_file)) // 4

        # rewind the handle to the beginning of the file.
        reads.seek(0)

        # process the reads (with progress!)
        with Progress() as progress:
            task = progress.add_task(
                "Loading reads", total=num_reads, unit="reads", width=40
            )
            for read in SeqIO.parse(reads, "fastq"):
                # print(f"read len: {len(read)}")
                kmers = create_possible_kmers(k, str(read.seq))
                # print(f"kmers len: {len(kmers)}")
                curr_graph = create_graph(curr_graph, kmers)

                revere_kmers = create_possible_kmers(
                    k, str(read.seq.reverse_complement())
                )
                # print(f"kmers len: {len(kmers)}")
                curr_graph = create_graph(curr_graph, revere_kmers)

                progress.update(task, advance=1)

    print(f"Processed {num_reads} reads.")

    return curr_graph


def main():
    print(create_possible_kmers(5, "AGCCGATCAT"))


if __name__ == "__main__":
    main()
