package gochujang

// borrowing heavily from Stephen Smith, but this is just for practice

type Node struct {
	label    string
	length   float64
	parent   *Node
	children []*Node
	data     map[string]int
	istip    bool
	note     string
	number   int
}

func NewNode() *Node {
	return &Node{}
}

// func parse_nwk(nwk string) (root *Node) {
// 	rt := Node{label: "", length: 0.0, parent: nil}
// 	for i, v := range nwk {
// 		if v == "(" {

// 		}
// 	}
// }

// ((A,B),C,D);

// borrowing heavily from https://talks.golang.org/2011/lex.slide#20
