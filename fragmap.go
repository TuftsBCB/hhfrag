package hhfrag

import (
	"fmt"
	"os"
	"runtime"
	"sort"
	"sync"

	"github.com/TuftsBCB/apps/hhsuite"
	"github.com/TuftsBCB/fragbag"
	"github.com/TuftsBCB/fragbag/bow"
	"github.com/TuftsBCB/io/fasta"
	"github.com/TuftsBCB/io/hmm"
	"github.com/TuftsBCB/seq"
)

type MapConfig struct {
	WindowMin       int
	WindowMax       int
	WindowIncrement int
	Blits           bool
}

var DefaultConfig = MapConfig{
	WindowMin:       30,
	WindowMax:       35,
	WindowIncrement: 5,
	Blits:           true,
}

func getOneFastaSequence(queryFasta string) (s seq.Sequence, err error) {
	fquery, err := os.Open(queryFasta)
	if err != nil {
		return
	}
	defer fquery.Close()

	seqs, err := fasta.NewReader(fquery).ReadAll()
	if err != nil {
		return
	} else if len(seqs) == 0 {
		err = fmt.Errorf("No sequences found in '%s'.", queryFasta)
		return
	} else if len(seqs) > 1 {
		err = fmt.Errorf("%d sequences found in '%s'. Expected only 1.",
			len(seqs), queryFasta)
		return
	}
	s = seqs[0]
	return
}

func (m MapConfig) MapFromFasta(pdbDb PDBDatabase, seqDb hhsuite.Database,
	queryFasta string) (*FragmentMap, error) {

	qseq, err := getOneFastaSequence(queryFasta)
	if err != nil {
		return nil, err
	}

	queryHHM, err := hhsuite.BuildHHM(
		hhsuite.HHBlitsDefault, hhsuite.HHMakePseudo, seqDb, queryFasta)
	if err != nil {
		return nil, err
	}
	return m.computeMap(pdbDb, qseq, queryHHM)
}

func (m MapConfig) MapFromHHM(pdbDb PDBDatabase, seqDb hhsuite.Database,
	queryFasta string, queryHHM string) (*FragmentMap, error) {

	qseq, err := getOneFastaSequence(queryFasta)
	if err != nil {
		return nil, err
	}

	fquery, err := os.Open(queryHHM)
	if err != nil {
		return nil, err
	}
	defer fquery.Close()

	qhhm, err := hmm.ReadHHM(fquery)
	if err != nil {
		return nil, err
	}
	return m.computeMap(pdbDb, qseq, qhhm)
}

func (m MapConfig) computeMap(
	pdbDb PDBDatabase, qseq seq.Sequence, qhhm *hmm.HHM) (*FragmentMap, error) {

	type maybeFrag struct {
		frags Fragments
		err   error
	}

	wg := new(sync.WaitGroup)
	jobs := make(chan int, 10)
	fragsChan := make(chan maybeFrag, 10)
	workers := runtime.GOMAXPROCS(0)
	if workers < 1 {
		workers = 1
	}

	for i := 0; i < workers; i++ {
		go func() {
			wg.Add(1)
			defer wg.Done()

			min, max := m.WindowMin, m.WindowMax
		CHANNEL:
			for start := range jobs {
				var best *Fragments
				for end := min; end <= max && (start+end) <= qseq.Len(); end++ {
					frags, err := FindFragments(
						pdbDb, m.Blits, qhhm, qseq, start, start+end)
					if err != nil {
						fragsChan <- maybeFrag{
							err: err,
						}
						continue CHANNEL
					}
					if best == nil || frags.better(*best) {
						best = frags
					}
				}
				fragsChan <- maybeFrag{
					frags: *best,
				}
			}
		}()
	}
	go func() {
		for s := 0; s <= qseq.Len()-m.WindowMin; s += m.WindowIncrement {
			jobs <- s
		}
		close(jobs)
		wg.Wait()
		close(fragsChan)
	}()

	fmap := &FragmentMap{
		Name:     qseq.Name,
		Segments: make([]Fragments, 0, 50),
	}
	for maybeFrag := range fragsChan {
		if maybeFrag.err != nil {
			return nil, maybeFrag.err
		}
		fmap.Segments = append(fmap.Segments, maybeFrag.frags)
	}
	sort.Sort(fmap)
	return fmap, nil
}

type FragmentMap struct {
	Name     string
	Segments []Fragments
}

func (fmap *FragmentMap) Len() int {
	return len(fmap.Segments)
}

func (fmap *FragmentMap) Less(i, j int) bool {
	return fmap.Segments[i].Start < fmap.Segments[j].Start
}

func (fmap *FragmentMap) Swap(i, j int) {
	fmap.Segments[i], fmap.Segments[j] = fmap.Segments[j], fmap.Segments[i]
}

func (fmap *FragmentMap) StructureBow(lib fragbag.StructureLibrary) bow.Bowed {
	bag := bow.NewBow(lib.Size())
	for _, fragGroup := range fmap.Segments {
		for _, frag := range fragGroup.Frags {
			bag = bag.Add(bow.StructureBow(lib, frag.CaAtoms))
		}
	}
	return bow.Bowed{Id: fmap.Name, Bow: bag}
}
