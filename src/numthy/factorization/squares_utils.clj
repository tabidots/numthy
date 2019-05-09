(ns numthy.factorization.squares-utils
  "Workhorse functions for the family of congruence-of-squares factorization algorithms."
  (:require [clojure.core.reducers :as r]
            [clojure.core.matrix :as m]
            [clojure.math.numeric-tower :refer [gcd expt abs]]
            [clojure.string :as s]
            [numthy.modular-arithmetic.utils :refer [mod-mul]]))

(defn smoothness-bound
  "Ceiling of e ^ (1/2 sqrt(ln n ln ln n))."
  ;; https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html
  ;; and some papers as well
  [n]
  (-> (*' (Math/log n) (Math/log (Math/log n)))
      Math/sqrt (* 1/2) Math/exp))

(defn smooth-factorization
  "Given an integer n, returns a map of the prime factors in primes that divide n
  and their exponents. Does not check if n factors completely over the primes."
  [n primes]
  (let [n' (abs n)]
    (->> (for [p primes]
           (loop [x n' exp 0]
             (cond
               (neg? p)         (if (neg? n) {-1 1} {-1 0})
               (pos? (mod x p)) {p exp}
               :else            (recur (/ x p) (inc exp)))))
         (apply merge))))

(defn factorize-rhs
  "Given a factor base and a hash-map of {x y} values where y is known to be
  smooth over the factor base, returns {x m} where m is an unsorted map of the
  prime factors of y and their exponents."
  [factor-base m]
  (reduce-kv (fn [relations x y]
               (assoc relations x (smooth-factorization y factor-base)))
             {} m))

(defn row->mod2
  "Given a vector, reduces all of its entries mod 2."
  [row]
  (mapv #(mod % 2) row))

(defn create-gf2-matrix
  "Given a factor base and collection of unsorted exponent maps, creates a GF(2) matrix
  (i.e., matrix where all entries are mod 2) from the maps."
  [base ms]
  (let [blank-map (into (sorted-map) (zipmap base (repeat 0)))]
    (mapv (comp row->mod2 vals #(merge blank-map %)) ms)))

(defn gf2-matrix->column-wise-binary
  "Encodes each column of a matrix in GF(2) as a binary number. Uses BigIntegers
  to accommodate matrices of arbitrary numbers of columns and rows."
  [m]
  (apply map (comp #(BigInteger. % 2)
                   #(s/reverse %)
                   #(reduce str %)
                   vector)
             m))

(comment
 "Keep as test data"
  (def m [[1 1 0 0]
          [1 1 0 1]
          [0 1 1 1]
          [0 0 1 0]
          [0 0 0 1]]))

(defn- lowest-set-bit
  "Finds the lowest set bit of a binary number. Only works for native Clojure numbers.
  This function is left here for reference."
  [a]
  (when (pos? a)
    (loop [b (.and a (- a)) low-bit -1]
      (if (zero? b) low-bit
        (recur (bit-shift-right b 1) (inc low-bit))))))

(defn get-pivot
  "Finds the pivot (first '1') in a binary-encoded column of a GF(2) matrix. Since
  the columns are encoded so that the top entry of the column is the rightmost bit,
  the first 1 is the lowest set bit of the number. When the pivot is -1, there is no
  pivot, so returns nil. Uses Java interop to handle BigIntegers."
  [n]
  (let [lb (.getLowestSetBit n)]
    (when-not (neg? lb)
      (long lb))))

(defn find-winning-subsets
  "Uses Koç & Arachchige's algorithm to perform fast Gaussian elimination in GF(2) and
  find all subsets of rows that add to the zero vector."
  ;; With reference to https://github.com/skollmann/PyFactorise/blob/master/factorise.py
  ;; Original paper: https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
  [m]
  (let [num-cols          (count (first m))
        num-rows          (count m)
        columns           (to-array (gf2-matrix->column-wise-binary m))
        ^booleans marked? (boolean-array num-rows false)
        ^longs    pivots  (long-array num-cols -1)]
    (dotimes [j num-cols]
             (when-some [i (get-pivot (aget columns j))]
               (aset-long pivots j i)
               (aset-boolean marked? i true)
               (dotimes [k num-cols]
                        (when (and (not= k j) (.testBit (aget columns k) i))
                          (aset columns k (.xor (aget columns k) (aget columns j)))))))
    (for [i (range num-rows) :when (not (aget ^booleans marked? i))]
      (cons i
            (for [j (range num-cols) :when (.testBit (aget columns j) i)]
              (aget ^longs pivots j))))))

(defn make-x-from
  "Given a list of indices, extracts the corresponding z-values from the candidates
  and multiplies them together to give x."
  [idxs zs n]
  (->> (m/order zs idxs)
       (r/fold (fn ([] 1) ([a b] (mod-mul a b n))))))

(defn make-y-from
  "Given a list of indices, extracts the corresponding exponent map and adds their values
  (which is equivalent to multiplying the integers they represent). This gives y^2.
  Then finds y by halving the entries of the sum vector, and converting that back to an integer."
  [idxs expm n]
  (->> (m/order expm idxs)
       (apply merge-with +)
       (r/fold (r/monoid (fn [a b] (mod-mul a b n)) (constantly 1))
               (fn [res b e]
                 (-> (expt b (quot e 2))
                     (mod-mul res n))))))

(defn find-factor
  "Given a hash-map of relations in the form {x m} where m is an exponent map of y,
  a transformed x-value, creates an exponent matrix over GF(2) from the maps of all ys,
  and solves for its left null space to find subsets of rows that add to a square. Then
  searches for a congruence of squares that yield a factor of n."
  [n factor-base relations]
  (let [xs            (vec (keys relations))
        exponent-maps (vec (vals relations))
        matrix        (create-gf2-matrix factor-base exponent-maps)]
    (->> (for [subset (find-winning-subsets matrix)]
           (let [x (make-x-from subset xs n)
                 y (make-y-from subset exponent-maps n)]
             (gcd (- x y) n)))
         (filter #(< 1 % n))
         (first))))
