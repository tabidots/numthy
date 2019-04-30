(ns numthy.factorization.squares-utils
  "Workhorse functions for the family of congruence-of-squares factorization algorithms."
  (:require [clojure.math.numeric-tower :refer [gcd expt]]
            [clojure.core.reducers :as r]
            [clojure.core.matrix :as m]
            [numthy.helpers :refer [divisible?]]
            [numthy.linear-algebra :as linalg]
            [numthy.modular-arithmetic.utils :refer [mod-mul]]))

(defn smoothness-bound
  "Ceiling of e ^ (1/2 sqrt(ln n ln ln n))."
  ;; https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html
  ;; and some papers as well
  [n]
  (-> (*' (Math/log n) (Math/log (Math/log n)))
      Math/sqrt (* 1/2) Math/exp))

(defn smooth?
  "Returns [exponent-vector] if n is smooth over the given empty hash-map of primes."
  [n primes prime-map]  ;; extra arg because iterating over keySeq is slow
  (loop [x      n
         m      prime-map
         primes primes]
    ;; TODO: Use an acc value to keep track of where we are in the primes vector.
    ;;       Implement Early Abort Strategy.
    ;;       https://books.google.com.vn/books?id=ITvaBwAAQBAJ&pg=PA202&lpg=PA202&dq=cfrac+multiplier+kn+%22continued+fraction%22+factorization+%22early+abort%22&source=bl&ots=6VvQXvUzSS&sig=ACfU3U0ldAojRlhZ_il97d_x9h_NF3Qp5g&hl=en&sa=X&ved=2ahUKEwj13N2ls_fhAhUQIIgKHWZ-B9EQ6AEwBXoECAkQAQ#v=onepage&q=cfrac%20multiplier%20kn%20%22continued%20fraction%22%20factorization%20%22early%20abort%22&f=false
    (if (= 1 x) (vec (vals m))
      (when-let [k (first primes)]
        (if (zero? (mod x k))
          (recur (/ x k) (update m k inc) primes)
          (recur x m (rest primes)))))))

(defn row->mod2
  "Given a vector, reduces all of its entries mod 2."
  [row]
  (mapv #(mod % 2) row))

(defn index-of-leading-1
  "Finds the index of the leading 1 in a row (vector). Returns nil if none."
  [row]
  (first (keep-indexed (fn [i x] (when (= 1 x) i)) row)))

(defn linearly-independent-idxs
  "Finds the subset of linearly independent rows in a matrix in ℚ."
  [m]
  (->> (m/transpose m)
       (linalg/reduced-row-echelon)
       (r/map index-of-leading-1)
       (r/remove nil?)
       (into [])))

(defn linearly-independent-idxs-gf2
  "Finds the subset of linearly independent rows in a matrix over GF(2)."
  [m]
  (->> (m/transpose m)
       (linalg/rref-gf2)
       (r/map index-of-leading-1)
       (r/remove nil?)
       (into [])))

(defn make-x-from-idxs
  "Given a list of indices, extracts the corresponding z-values from the candidates
  and multiplies them together to give x."
  [idxs zs n]
  (->> (m/order zs idxs)
       (r/fold (fn ([] 1) ([a b] (mod-mul a b n))))))

(defn make-y-from-idxs
  "Given a list of indices, extracts the corresponding rows from the original
  exponent matrix and adds their values (which is equivalent to multiplying the
  integers they represent). This gives y^2. Then finds y by halving the entries
  of the sum vector, and converting that back to an integer."
  [idxs expm factor-base n]
  (->> (m/order expm idxs)
       (apply m/add)
       (zipmap factor-base)
       (r/fold (r/monoid (fn [a b] (mod-mul a b n)) (constantly 1))
               (fn ([res b e]
                    (-> (expt b (quot e 2))
                        (mod-mul res n)))))))

(defn find-factor-from-congruent-relations
  "Separates a matrix of relations into linearly independent and dependent rows.
  For the submatrix of linearly independent rows A, and each dependent row b,
  solves the matrix equation Ax = b in GF(2) to find the vectors of the null space of A
  until a congruence of squares is found."
  [n factor-base relations]
  (let [zs               (vec (keys relations))
        expm             (vec (vals relations))
        binary-expm      (mapv row->mod2 expm)
        independent-idxs (linearly-independent-idxs-gf2 binary-expm)
        At               (-> binary-expm
                             (m/order 0 (vec independent-idxs))
                             (m/transpose))]
    (loop [dependent-idxs (->> (range (count expm))
                               (r/remove (set independent-idxs))
                               (r/foldcat))]
      (if-let [d-idx (first dependent-idxs)]
        (let [b    (get binary-expm d-idx)
              ;; Solve xA = b, that is A^T • x = b,  for each dependent row b.
              idxs (->> (linalg/solve-gf2 At b) ; old: (row->mod2 (linalg/least-squares At b))
                        ;; Example: If the independent indxs are [1 2 3 4 5] and the solution is [0 1 1 0 1],
                        ;; then congruent indxs are [- 2 3 - 5] -> [2 3 5].
                        (keep-indexed (fn [i x] (when-not (zero? x)
                                                  (get independent-idxs i))))
                        (cons d-idx)   ;; Include the index for row b. e.g. (10 2 3 5) if solving for row 10
                        (vec))
              x    (make-x-from-idxs idxs zs n)
              y    (make-y-from-idxs idxs expm factor-base n)
              fact (gcd (+ x y) n)]
          (if (and (not= (mod x n) (mod y n)) (< 1 fact n)) fact
            (recur (rest dependent-idxs))))
        "No factors found :("))))
