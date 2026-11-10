import argparse
from treelib import Tree, exceptions

def get_log_line_info(line):
	if not line.startswith("   "):
		raise ValueError("Expected a backtrace line, but got", line)
	infos = line.split(" ")[3:]
	last_bt = infos[0] == "at"
	addr = infos[1][:-1]
	func = infos[2]
	debug = infos[3][1:-1] if infos[3] != "(in" else infos[4][:-1]
	return last_bt, addr, func, debug


def get_node_infos(line, parent_id):
	if not parent_id :
		parent_id = ""
	if parent_id.startswith("0x"):
		parent_id = parent_id[2:]
	infos = get_log_line_info(line)
	node_id = hex(hash(int(infos[1] + parent_id, 16)))

	return node_id, infos


def build_errors_tree(backtraces, depth=6):
	# Call tree keeping track of the errors. Each nodes as the following
	# attributes:
	#	* Number of errors in the sub tree
	tree = Tree()

	for bt in backtraces:
		bt.reverse()
		parent_id = None
		for i in range(min(len(bt), depth)):
			try:
				node_id, infos = get_node_infos(bt[i], parent_id)
			except ValueError:
				continue
			node = tree.get_node(node_id)
			if not node:
				label = "in {} ({})".format(infos[2], infos[3])
				data = 0 if infos[0] else 1
				try:
					tree.create_node(label, node_id, parent_id, data=data)
				except exceptions.MultipleRootError:
					print("MultipleRootError(A tree takes one root merely.)")
					print(bt[i])
					print(bt[i])
					print(label, node_id, parent_id, data)
					print(tree)
					exit(-1)
			else:
				node.data += 1
			parent_id = node_id
	return tree


def print_backend(tree, nid=None, level=0, idhidden=True, filter=None,
					key=None, reverse=False, line_type='ascii-ex',
					data_property=None, func=print):
	"""
	Another implementation of printing tree using Stack
	Print tree structure in hierarchy style.
	For example:
	.. code-block:: bash
		Root
		|___ C01
		|    |___ C11
		|         |___ C111
		|         |___ C112
		|___ C02
		|___ C03
		|    |___ C31
	A more elegant way to achieve this function using Stack
	structure, for constructing the Nodes Stack push and pop nodes
	with additional level info.
	UPDATE: the @key @reverse is present to sort node at each
	level.
	"""
	# Factory for proper get_label() function
	def get_label(node):
		if data_property:
			if idhidden:
				return getattr(node.data, data_property)
			else:
				return "%s[%s]" % (getattr(node.data, data_property), node.identifier)
		else:
			if idhidden:
				return node.tag
			else:
				return "%s[%s]" % (node.tag, node.identifier)

	# legacy ordering
	if key is None:
		def key(node):
			return node

	# iter with func
	for pre, node in __get(tree, nid, level, filter, key, reverse,
								line_type):
		label = get_label(node)
		nerrs = ' + ' if node.is_leaf() and node.data == 0 else str(node.data).center(3, ' ')
		func('{0}{1}'.format(pre, label).encode('utf-8'), nerrs)

def __get(tree, nid, level, filter_, key, reverse, line_type):
	# default filter
	if filter_ is None:
		def filter_(node):
			return True

	# render characters
	dt = {
		'noguides': (' ', '    ', '    '),
		'ascii': ('|', '|-- ', '+-- '),
		'ascii-ex': ('\u2502', '\u251c\u2500\u2500 ', '\u2514\u2500\u2500 '),
		'ascii-exr': ('\u2502', '\u251c\u2500\u2500 ', '\u2570\u2500\u2500 '),
		'ascii-em': ('\u2551', '\u2560\u2550\u2550 ', '\u255a\u2550\u2550 '),
		'ascii-emv': ('\u2551', '\u255f\u2500\u2500 ', '\u2559\u2500\u2500 '),
		'ascii-emh': ('\u2502', '\u255e\u2550\u2550 ', '\u2558\u2550\u2550 '),
	}[line_type]

	return __get_iter(tree, nid, level, filter_, key, reverse, dt, [])

def __get_iter(tree, nid, level, filter_, key, reverse, dt, is_last):
	dt_vertical_line, dt_line_box, dt_line_corner = dt

	nid = tree.root if nid is None else nid
	node = tree[nid]

	if level == 0:
		yield "", node
	else:
		leading = ''.join(map(lambda x: dt_vertical_line + ' ' * 3
		if not x else ' ' * 2, is_last[0:-1]))
		lasting = dt_line_corner if is_last[-1] else dt_line_box
		yield leading + lasting, node

	if filter_(node) and node.expanded:
		children = [tree[i] for i in node.successors(tree._identifier) if filter_(tree[i])]
		idxlast = len(children) - 1
		if key:
			children.sort(key=key, reverse=reverse)
		elif reverse:
			children = reversed(children)
		level += 1
		for idx, child in enumerate(children):
			is_last.append(idx == idxlast)
			for item in __get_iter(tree, child.identifier, level, filter_,
										key, reverse, dt, is_last):
				yield item
			is_last.pop()

def show(tree, nid=None, level=0, idhidden=True, filter=None,
			key=None, reverse=False, line_type='ascii-ex'):
	tree.reader = ""

	def write(line, n_err='+', len_line=90):
		trailing_spaces = ' ' * (len_line-len(line.decode('utf-8'))-len(n_err))
		tree.reader += line.decode('utf-8') + trailing_spaces + n_err + "\n"

	print_backend(tree, nid, level, idhidden, filter,
								key, reverse, line_type, func=write)


parser = argparse.ArgumentParser(description="Pretty display of valgrind's error reports")
parser.add_argument('valgrind_reports', metavar='filename', type=str, nargs='+',
                    help="log file containing valgrind's report")
parser.add_argument("--outdir", "-o", type=str, default=None, help="path where to save the final reports")
parser.add_argument("--depth", "-d", type=int, default=6, help="depth of the call tree")
args = parser.parse_args()

# Process every report file given as input, poducing a call tree for each one
trees = []
for log_file in args.valgrind_reports:
	with open(log_file) as fp:
		report = fp.read().split("\n\n")
	backtraces = [bt.split("\n") for bt in report[1:]]
	print(log_file)
	t = build_errors_tree(backtraces, depth=args.depth)
	trees.append(t)

# Now we have all tree, we build a tree containing all common occurences, and different trees for all subset of trees containing a specific occurence
all_tree = Tree()
common_tree = Tree()
diff_tree = {}
for i in range(len(trees)):
	tree = trees[i]
	for node in tree.all_nodes():
		if all_tree.contains(node.identifier):
			continue
		all_tree.add_node(node, tree.parent(node.identifier))
		tree_containing_n = []
		for j in range(len(trees)):
			if trees[j].contains(node.identifier):
				tree_containing_n.append(args.valgrind_reports[j])
		if len(tree_containing_n) == len(trees):
			common_tree.add_node(node, tree.parent(node.identifier))
		else:
			tree_id = '-'.join(sorted([x.split("/")[-1].split("_")[-1].split('.')[0] for x in tree_containing_n]))
			if tree_id not in diff_tree:
				diff_tree[tree_id] = Tree()
			nodes_to_add = [node]
			tmp_node = node
			while not diff_tree[tree_id].contains(tree.parent(tmp_node.identifier)) and not tmp_node.is_root():
				nodes_to_add = [tree.parent(tmp_node.identifier)] + nodes_to_add
				tmp_node = tree.parent(tmp_node.identifier)
			for k in range(len(nodes_to_add)):
				try:
					diff_tree[tree_id].add_node(nodes_to_add[k], tree.parent(nodes_to_add[k].identifier))
				except:
					continue

try:
	show(all_tree, line_type=u"noguides")
	if args.outdir:
		output_filename = f"{args.outdir}_{args.depth}-all_errors.log"
		with open(output_filename, 'w') as fp:
			fp.write(all_tree.reader)
	else:
		print(common_tree.reader)
except Exception as e:
	print(f"No errors ({e})")

if len(diff_tree) > 0:
	try:
		show(common_tree, line_type=u"noguides")
		if args.outdir:
			output_filename = f"{args.outdir}/{args.depth}-common_errors.log"
			with open(output_filename, 'w') as fp:
				fp.write(common_tree.reader)
		else:
			print(common_tree.reader)
	except Exception:
		print("No common errors")

	for t in diff_tree:
		show(diff_tree[t], line_type=u"noguides")
		if args.outdir:
			output_filename = f"{args.outdir}/{args.depth}-subset-{t}.log"
			with open(output_filename, 'w') as fp:
				fp.write(diff_tree[t].reader)
		else:
			print(diff_tree[t].reader)
