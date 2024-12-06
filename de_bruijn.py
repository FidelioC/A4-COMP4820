import graphviz
import sys, time
from collections import deque
from rich.progress import Progress
from Bio import SeqIO
from Bio.Seq import Seq
import os
import click

from resource import getrusage, RUSAGE_SELF

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
    # if kmer1 exist, we want to make a new node (kmer2), and points kmer1 to kmer2
    if kmer_2 not in graph[kmer_1]:
        graph[kmer_1][kmer_2] = 1
    else:
        graph[kmer_1][kmer_2] += 1
    graph[kmer_2] = {}


def edit_graph_kmer2_exist(graph, kmer_1, kmer_2):
    # if kmer2 exist, we want to make a new node, and kmer1 points to kmer2
    graph[kmer_1] = {}
    graph[kmer_1][kmer_2] = 1


def edit_graph_both_dont_exist(graph, kmer_1, kmer_2):
    graph[kmer_1] = {}
    graph[kmer_1][kmer_2] = 1
    graph[kmer_2] = {}


def create_graph(curr_graph, all_kmers, all_kmers_reverse):
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
    count = 0
    while index + 1 < len(all_kmers):
        if count % 2 == 0:
            kmer_1 = all_kmers[index]
            kmer_2 = all_kmers[index + 1]
        else:
            kmer_1 = all_kmers_reverse[index]
            kmer_2 = all_kmers_reverse[index + 1]
            index += 1

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
        count += 1
        # print(f"kmer1: {kmer_1}, kmer2: {kmer_2}, graph: {graph}\n")

        # dot.edge(kmer_1, kmer_2, label=str(weight))
        # print(graph)

    return graph


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
        curr_parent_outdegree = get_node_outdegree(graph, node)
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
def remove_bubble(graph, k):
    total_bubbles_removed = 0
    # identify bubbles w/ bfs
    potential_bubbles = find_bubble(graph)
    # loop until there's no more bubbles
    while len(potential_bubbles) > 0:
        # print(f"potential bubbles: {potential_bubbles}")
        for bubble in potential_bubbles:
            # find edge w/ lowest weight
            min_child = find_lowest_weight_children(graph[bubble])
            # print(f"bubble: {bubble}, min_child: {min_child}")

            # delete the edge w/ lowest weight
            graph = remove_node_from_graph(min_child, graph)
            # print(f"removing bubble at: {min_child}")
            total_bubbles_removed += 1

            # prune any tips
            graph, total_prune_outdegree, total_prune_indegree = remove_graph_all_tips(
                graph, k
            )
            # print("remove_graph_all_tips finished")
            total_bubbles_removed += total_prune_indegree + total_prune_outdegree
        # print(f"graph after removing bubble: {graph}")
        # repeat until no bubbles exist
        potential_bubbles = find_bubble(graph)
        # print("remove potential bubbles finished")
    # print("remove bubbles finished")

    return graph, total_bubbles_removed


def find_lowest_weight_children(children: dict):
    """
    return the children with the lowes weight
    output: e.g., "ATCAT"
    """
    min_weight = sys.maxsize
    min_node = None
    curr_weight = 0
    for children, weight in children.items():
        curr_weight = weight
        if curr_weight < min_weight:
            min_weight = curr_weight
            min_node = children

    return min_node


def find_bubble(graph):
    """
    find potential bubbles in nodes with outdegree > 1
    return value: dict of bubbles, where key: (startbubble node, end bubble node)
    """
    potential_bubbles = {}
    result = None
    outdegree_greater_zero = find_all_outdegree_greater_zero(graph)
    # print(f"outdegree_greater_zero {outdegree_greater_zero}")
    for node in outdegree_greater_zero:
        if node in graph:
            # print(f"\nNode: {node}\n")
            result = traverse_find_bubble_iterative(graph, node)
        if result:
            potential_bubbles[node] = result

        result = None

    return potential_bubbles


def traverse_find_bubble_iterative(graph, start_node):
    visited = set()
    queue = deque([start_node])  # Initialize queue with start node
    end_node = None

    while queue:  # until queue is empty
        # print(queue)
        node = queue.popleft()  # next node

        if node in visited:  # node has been revisited
            if node == start_node:
                # if we cycle back to the start_node, it's not a valid bubble
                return None
            end_node = node
            return (start_node, end_node)  # found bubble

        visited.add(node)  # Mark the node as visited

        if get_node_outdegree(graph, node) > 1 and node != start_node:
            return None

        # enqueue unvisited neighbors
        for neighbor in graph[node]:
            if neighbor not in visited:
                queue.append(neighbor)

    return None  # no bubble found


def find_all_outdegree_greater_zero(graph):
    """
    find all nodes that have outdegree greater than one
    """
    outdegree_greater_zero = []
    for node, neighbours in graph.items():
        if len(neighbours) > 1:
            outdegree_greater_zero.append(node)

    return outdegree_greater_zero


def get_node_outdegree(graph, node):
    return len(graph[node])


def remove_node_from_graph(node, graph):
    """given a node, remove it from the graph,
    first will delete edges that's pointing to this node,
    then, delete the node itself
    """
    parents = find_parent(graph, node)
    for parent in parents:
        del graph[parent][node]  # delete the parent's edge pointing to the node
    del graph[node]  # delete the node itself

    return graph


# =============== all errors ================
def remove_all_errors(graph, k):
    graph, total_prune_outdegree, total_prune_indegree = remove_graph_all_tips(graph, k)
    graph, total_bubbles_removed = remove_bubble(graph, k)

    return graph, (total_prune_indegree + total_prune_outdegree), total_bubbles_removed


# =============== contigs ===================
def create_contigs(graph: dict):
    contigs = {}
    contig_id = 1
    indegree_zero = find_all_indegree_zero(graph)

    for node in indegree_zero:
        curr_node = node
        contig = curr_node
        visited_in_contig = set()

        # traverse nodes and also handle cycle
        while (
            get_node_outdegree(graph, curr_node) >= 1
            and curr_node not in visited_in_contig
        ):
            visited_in_contig.add(curr_node)
            neighbours = graph[curr_node]

            # here, if there's still a node with outdegree > 1,
            # I will arbitrarily pick the highest weight.
            next_node = max(neighbours, key=neighbours.get)
            contig += next_node[-1]  # add last char to contig
            curr_node = next_node

        contigs[contig_id] = contig

        contig_id += 1

    return contigs


def write_contigs(contigs: dict):
    for id, contig in contigs.items():
        filename = os.path.join(f"{id}-covid.fq")
        with open(filename, "w") as fastq_file:
            fastq_file.write(f"@{id}-covid\n")
            fastq_file.write(f"{contig}\n")  # sequence
            fastq_file.write("+\n")


# =============== file/ read processing ==================
def create_graph_list_reads(all_reads, k):
    curr_graph = {}
    # dot = graphviz.Digraph(format="svg")
    for read in all_reads:
        # print(f"============= CURRENT READ {read} ============= \n")
        # print(f"graph: {curr_graph}")
        curr_graph = create_graph_two_read(curr_graph, read, k)
        # print(f"KMERS: {kmers}")
    # dot.render("de_bruijn_graph_output")
    return curr_graph


def get_reverse_complement_kmers(k, sequence):
    rev_comp_sequence = str(Seq(sequence).reverse_complement())
    return create_possible_kmers(k, rev_comp_sequence)


def create_graph_two_read(graph, read, k):
    to_seq = Seq(read)

    kmers = create_possible_kmers(k, str(to_seq))
    # print(f"kmers len: {len(kmers)}")
    # graph = create_graph(graph, kmers)

    reverse_kmers = get_reverse_complement_kmers(k, str(to_seq))
    # print(f"kmers len: {len(kmers)}")
    graph = create_graph(graph, kmers, reverse_kmers)

    return graph


def create_graph_file_reads(read_file, k):
    curr_graph = {}
    len_read = 0
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
                len_read += len(read)
                create_graph_two_read(curr_graph, read.seq, k)
                progress.update(task, advance=1)

            # print(f"read len: {len_read}")

    # print(f"Processed {num_reads} reads.")

    return curr_graph, num_reads


# =============== print output ===========================
def mem_usage():
    # 'maxrss' is the maximum resident set size
    max_memory = getrusage(RUSAGE_SELF).ru_maxrss
    print(f"* Maximum resident set size (memory usage): {max_memory}kb")


def print_output(
    total_read,
    init_nodes,
    init_input_nodes,
    init_output_nodes,
    pruned_tips,
    input_nodes_remain,
    output_nodes_remain,
    removed_bubbles,
    contigs: dict,
):
    print("Output")
    print("======")

    print(f"\n* Added {total_read} reads to the graph.")
    print(f"* Graph has:")
    print(f"\t* {init_nodes} total nodes.")
    print(f"\t* {init_input_nodes} input nodes.")
    print(f"\t* {init_output_nodes} output nodes.")

    print(f"\n* Pruned {pruned_tips} nodes.")
    print(f"\t* {input_nodes_remain} input nodes remain.")
    print(f"\t* {output_nodes_remain} output nodes remain.")
    print(f"\t* Removed {removed_bubbles} nodes in bubbles.")
    print(f"* Reconstructed {len(contigs)} sequences:")

    for id, contig in contigs.items():
        print(f"* {id}-covid.fq ({len(contig)}bp)")

    mem_usage()


def print_settings(reads, kmer_size):
    print("Settings")
    print("========")

    print(f"\n* Reads: {reads}")
    print(f"* k:     {kmer_size}\n")


@click.command()
@click.option("--reads", required=True)
@click.option("--kmer-size", default=37)
def main(reads, kmer_size):
    k_file = kmer_size
    file_name = reads

    print_settings(reads, kmer_size)

    graph, total_read = create_graph_file_reads(file_name, k_file)

    init_nodes = len(graph)
    init_input_nodes = len(find_all_indegree_zero(graph))
    init_output_nodes = len(find_all_outdegree_zero(graph))
    # visualize_graph(graph, "./test_graph_main_beforeerror")
    graph, pruned_tips, removed_bubbles = remove_all_errors(graph, k_file)
    # print("remove all errors finished")
    # visualize_graph(graph, "./test_graph_main")

    output_nodes_remain = len(find_all_outdegree_zero(graph))
    input_nodes_remain = len(find_all_indegree_zero(graph))

    contigs = create_contigs(graph)

    write_contigs(contigs)

    print_output(
        total_read,
        init_nodes,
        init_input_nodes,
        init_output_nodes,
        pruned_tips,
        input_nodes_remain,
        output_nodes_remain,
        removed_bubbles,
        contigs,
    )


if __name__ == "__main__":
    main()
