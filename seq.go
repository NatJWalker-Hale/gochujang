package gochujang

import (
	"bufio"
	"log"
	"os"
	"strings"
)

type DataType string

// datatype constants
const (
	Nucleotide DataType = "nuc"
	AminoAcid  DataType = "aa"
	MultiState DataType = "mult"
)

type Sequence struct {
	alphabet DataType
	name     string
	sequence string
	BF       []float64
	gc       float64
}

type SequenceDB struct {
	alphabet  DataType
	sequences []*Sequence
	aligned   bool
	length    int
	BF        []float64
}

func GetStates(alphabet DataType) []string { // helper function for state constants
	if alphabet == "nuc" {
		return []string{"A", "T", "G", "C"}
	} else if alphabet == "aa" {
		return []string{"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"}
	} else {
		return nil // placehold here for multistate
	}
}

func NewSequence() *Sequence {
	return &Sequence{}
}

func NewSequenceDB() *SequenceDB {
	return &SequenceDB{}
}

func ReadSeqsFromFile(path string) (seqs SequenceDB) {
	file, err := os.Open(path)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	first := true
	var cname string
	var cseq string
	for scanner.Scan() {
		if scanner.Text()[0:1] == ">" {
			if first {
				first = false
				cname = scanner.Text()[1:] // read first name
			} else {
				seq := NewSequence() // yield last entry
				seq.name = cname
				seq.sequence = cseq
				seqs.sequences = append(seqs.sequences, seq)
				cseq = ""
				cname = scanner.Text()[1:] // read new name
			}
		} else {
			cseq += scanner.Text() // concat multiple lines if present
		}
	}
	seq := NewSequence() // get last entry
	seq.name = cname
	seq.sequence = cseq
	seqs.sequences = append(seqs.sequences, seq)

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	for _, s := range seqs.sequences {
		s.GuessAlphabet()
		s.CalcBF()
	}
	alph := seqs.sequences[0].alphabet
	seqs.alphabet = alph
	for _, s := range seqs.sequences {
		if s.alphabet != alph {
			log.Fatal("sequences are not of the same alphabet!")
			// all seqs in DB should be of same alphabet
		}
	}
	seqs.aligned = true
	seqlen := len(seqs.sequences[0].sequence)
	for _, s := range seqs.sequences {
		if len(s.sequence) != seqlen {
			seqs.aligned = false
		}
	}
	if seqs.aligned {
		seqs.length = seqlen
	}
	seqs.CalcBF()
	return
}

func (s Sequence) GetFasta() string {
	return ">" + s.name + "\n" + s.sequence
}

func (s SequenceDB) GetFasta() (out string) {
	for _, v := range s.sequences {
		out += v.GetFasta() + "\n"
	}
	return
}

func (s *Sequence) GuessAlphabet() {
	dna := map[string]int{
		"A": 0,
		"T": 1,
		"G": 2,
		"C": 3,
		"-": 4,
		"N": 4,
	} // add extended IUPAC later
	s.alphabet = "nuc"
	for _, v := range s.sequence {
		if _, exists := dna[string(v)]; exists {
			continue
		} else {
			s.alphabet = "aa"
		}
	}
}

func (s *Sequence) CalcBF() {
	if s.alphabet == "nuc" {
		NUCs := GetStates(s.alphabet)
		NUCcount := make(map[string]int)
		NUCprop := make(map[string]float64)
		tot := 0
		for _, n := range NUCs {
			NUCcount[n] = strings.Count(s.sequence, n)
			tot += NUCcount[n]
		}
		for _, n := range NUCs {
			NUCprop[n] = float64(NUCcount[n]) / float64(tot)
			s.BF = append(s.BF, NUCprop[n])
		}
		s.gc = s.BF[2] + s.BF[3]
	} else if s.alphabet == "aa" {
		AAs := GetStates(s.alphabet)
		AAcount := make(map[string]int)
		AAprop := make(map[string]float64)
		tot := 0
		for _, a := range AAs {
			AAcount[a] = strings.Count(s.sequence, a)
			tot += AAcount[a]
		}
		for _, a := range AAs {
			AAprop[a] = float64(AAcount[a]) / float64(tot)
			s.BF = append(s.BF, AAprop[a])
		}
	} else {
		s.GuessAlphabet()
		s.CalcBF()
	}
}

func (s *SequenceDB) CalcBF() {
	if s.alphabet == "nuc" {
		NUCs := GetStates(s.alphabet)
		NUCcount := make(map[string]int)
		NUCprop := make(map[string]float64)
		tot := 0
		for _, v := range s.sequences {
			for _, n := range NUCs {
				NUCcount[n] += strings.Count(v.sequence, n)
				tot += strings.Count(v.sequence, n) // don't reinclude previous sequences counts
			}
		}
		for _, n := range NUCs {
			NUCprop[n] = float64(NUCcount[n]) / float64(tot)
			s.BF = append(s.BF, NUCprop[n])
		}
	} else {
		AAs := GetStates(s.alphabet)
		AAcount := make(map[string]int)
		AAprop := make(map[string]float64)
		tot := 0
		for _, v := range s.sequences {
			for _, a := range AAs {
				AAcount[a] += strings.Count(v.sequence, a)
				tot += strings.Count(v.sequence, a)
			}
		}
		for _, a := range AAs {
			AAprop[a] = float64(AAcount[a]) / float64(tot)
			s.BF = append(s.BF, AAprop[a])
		}
	}
}

func (s SequenceDB) GetColumns() map[int]string { // iterate through and populate map of position, column string
	if !s.aligned {
		log.Fatal("cannot return columns, sequences are not aligned!")
	}
	columns := make(map[int]string)
	for pos := 0; pos < s.length; pos++ {
		for _, v := range s.sequences {
			columns[pos] += string(v.sequence[pos])
		}
	}
	return columns
}

// func main() {
// 	aln := flag.String("s", "", "your sequences, in FASTA")
// 	flag.Parse()

// 	seqs := ReadSeqsFromFile(*aln)

// 	fmt.Println(seqs.aligned)
// 	fmt.Println(seqs.alphabet)

// 	fmt.Printf("%.4f\n", seqs.BF)

// 	fmt.Println(seqs.GetFasta())

// 	columns := seqs.GetColumns()
// 	fmt.Println(len(columns))
// 	for k, v := range columns {
// 		fmt.Println(k, v)
// 	}
// }
