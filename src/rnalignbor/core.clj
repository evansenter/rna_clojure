(ns rnalignbor.core)

(defn reload_env
  "Reload the environment"
  []
  (use `rnalignbor.core :reload-all))

(def temperature (+ 37 273.15))
(def boltzmann_constant 0.0019872370936902486) ; kcal / mol / K
(def base_pair_energy Math/exp (/ 1 (* boltzmann_constant temperature)))
(def min_loop_size 3)

(defn rand_sequence 
	"Generate a random RNA sequence of specified length"
	[length]
	(repeatedly length #(rand-nth [\a \u \g \c])))

(defn pairable
  "Boolean to determine if two characters form a valid base pair"
  [a b]
  (boolean (some #{b} ({ \a [\u] \u [\a \g] \g [\c \u] \c [\g] } a))))

(defn pairable_in_seq
  "Boolean to determine if two indices form a valid base pair"
  [rna_seq i j]
  (and 
    (> (- j i) min_loop_size)
    (pairable (nth rna_seq i) (nth rna_seq j))))

(defn _num_structures 
  "Calculates the number of structures for a given RNA sequence (doesn't include base triplets or pseudoknots)"
  [rna_seq & [{ :keys [i j] :or { i 0 j (-> rna_seq count dec) } }]] 
  (if (< (dec (- j i)) min_loop_size)
    1
    (+ 
      (_num_structures rna_seq {:i i :j (dec j) }) 
      (if (pairable_in_seq rna_seq i j) (_num_structures rna_seq { :i (inc i) :j (dec j) }) 0)
      (reduce 
        (fn [sum k] (+ sum (* 
          (_num_structures rna_seq { :i i :j (dec k) }) 
          (_num_structures rna_seq { :i (inc k) :j (dec j) }))))
        0
        (filter #(pairable_in_seq rna_seq % j) (range (inc i) (- j min_loop_size)))))))

(def num_structures (memoize _num_structures))