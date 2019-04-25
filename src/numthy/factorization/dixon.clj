(ns numthy.factorization.dixon
  (:require [clojure.math.numeric-tower :refer [gcd expt]]
            [numthy.helpers :refer [divisible?]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.modular-arithmetic.utils :refer [mod-pow mod-mul]]
            [numthy.perfect-powers :refer [perfect-power?]]
            [numthy.linear-algebra :as linalg]))

(defn- smoothness-bound
  "Ceiling of e ^ (1/2 sqrt(ln n ln ln n))."
  [n]
  (-> (*' (Math/log n) (Math/log (Math/log n)))
      Math/sqrt (* 1/2) Math/exp))

(defn- smooth?
  "Returns [exponent-vector] if n is smooth over the given list of primes."
  [primes n]
  (loop [x         n
         primes    primes
         exponents (->> (map #(hash-map % 0) primes)
                        (reduce into (sorted-map)))]
    (when-let [k (first primes)] ;; "nils out" if there are no more primes; i.e., n is not b-smooth
      (cond
        (= 1 x)            (vec (vals exponents))
        (divisible? x k) (recur (/ x k) primes (update exponents k inc))
        :else              (recur x (rest primes) exponents)))))

(defn- square?
  "Given an exponent vector, returns true if all its members are even."
  [expv]
  (every? even? expv))

(defn- row->mod2
  "Given a vector, reduces all of its entries mod 2."
  [row]
  (mapv #(mod % 2) row))

(defn- relations
  "Given a factor base (vector of primes) with B elements, accumulate B+1 values of z,
  between sqrt(n) and n, s.t. z^2 mod n is smooth over the factor base. Returns a hash-map
  with the zs as a vector and the exponent vectors in a vector (effectively an exponent matrix).
  Returns nil if no relations found."
  [n factor-base]
  (let [r (->> (drop (Math/sqrt n) (range (inc n)))
               (keep (fn [z]
                       (when-let [expv (smooth? factor-base (mod-pow z 2 n))]
                         {:z z :expv expv})))
               (take (inc (count factor-base))))]
    (when (not-empty r)
      {:zs (mapv :z r) :exponent-matrix (mapv :expv r)})))

(defn- segregate-rows
  "Segregates the rows of an exponent matrix, reduced mod 2, into:
  - A hash-map of its linearly independent rows (in the form {idx [exp-vec]})
  - A hash-map of its linearly dependent rows (in the form {idx [exp-vec]})
  This redundancy is due to the fact that the original, non mod-reduced matrix
  is also needed to find the congruences."
  [matrix]
  (let [expm->mod2          (mapv row->mod2 matrix)
        ;; ↓ "To find the linearly independent rows of the exponent matrix,"
        ;; "we can get the RREF of the exponent matrix transposed"
        crazy-matrix        (->> expm->mod2
                                 linalg/transpose
                                 linalg/reduced-row-echelon)
        ;; ↓ The idxs of leading 1s tell us the idxs of linearly independent
        ;; rows in the exponent matrix.
        leading-1-idx       (fn [row]
                              (first (keep-indexed (fn [i x]
                                                     (when (= 1 x) i))
                                                   row)))
        independent-idxs    (keep leading-1-idx crazy-matrix)
        dependent-idxs      (remove (set independent-idxs)
                                    (range (count matrix)))]
    {:independent-rows      (reduce (fn [res i]
                                      (assoc res i (get expm->mod2 i)))
                                    {} independent-idxs)
     :dependent-rows        (reduce (fn [res i]
                                      (assoc res i (get expm->mod2 i)))
                                    {} dependent-idxs)}))

(defn- find-congruent-sets
  "Given a matrix segregated into linearly independent and dependent rows,
  where A is the matrix of linearly independent rows, and b is a dependent row,
  this solves Ax = b for each b, mod 2. x is a vector containing 0s and 1s,
  where 1s mean keep the corresponding row in A, and 0s mean to drop it.
  Returns the idxs of the filtered rows in A with idx of b appended."
  [segregated-matrix]
  (let [A       (vec (vals (:independent-rows segregated-matrix)))
        indeps  (keys (:independent-rows segregated-matrix))
        ;; Winners are a solution to the matrix equation Ax = b and look like [0 1 1 0 1]
        ;; Returns solution or nil if all zeros.
        winners (fn [b]
                  (let [soln (->> (linalg/least-squares (linalg/transpose A) b)
                                  (row->mod2))]
                    (when (not-every? zero? soln) soln)))]
    (->> (:dependent-rows segregated-matrix)
         (reduce-kv (fn [res idx b]  ;; Solve Ax = b for each dependent row b.
                      (if-let [w (winners b)]
                        ;; If the independent indxs are [1 2 3 4 5] and winners are [0 1 1 0 1],
                        ;; then congruent indxs are [- 2 3 - 5] -> [2 3 5].
                        (->> (map (fn [a b] (when-not (zero? b) a)) indeps w)
                             (remove nil?)
                             (cons idx)  ;; Add the index for row b. e.g. (10 2 3 5) if solving for row 10
                             (conj res)) ;; Append this list of idxs to result.
                        res)) []))))

(defn- find-factors
  "Given x, y, and modulus n, uses congruence of squares mod n to determine whether
  gcd(x + y, n) and gcd(x - y, n) are factors of n. This function does not exclude
  trivial factors (1 and n)."
  [x y n]
  (when (and (not= (mod x n) (mod y n))
             (= (mod (*' x x) n) (mod (*' y y) n)))
    [(gcd (+ x y) n) (gcd (- x y) n)]))

(defn good-candidate?
  ;; http://micsymposium.org/mics_2011_proceedings/mics2011_submission_28.pdf
  [n]
  (when-not (or (.isProbablePrime (biginteger n) 5)
                (->> (range 2 (Math/log10 n))
                     (filter prime?)
                     (some #(divisible? n %)))
                (perfect-power? n))
    n))

(defn dixon-factorize
  "Uses Dixon's factorization algorithm to factorize a large integer n. Does not
  necessarily factorize n completely. Returns nil if no relations are found.
  Speed of algorithm depends on ability to find relations quickly. This is unpredictable,
  but generally does not happen."
  ;; With reference to https://crypto.stanford.edu/cs359c/17sp/projects/BrendonGo.pdf
  [n] ; try 16850989; 1078766140548460
  (let [B               (smoothness-bound n)
        factor-base     (filter prime? (range B))]
    (when-let [{zs :zs exponent-matrix :exponent-matrix} (relations n factor-base)]
      (let [congruent-sets  (->> exponent-matrix segregate-rows find-congruent-sets)
            x-y-pair        (fn [c-set]
                              ;; Given a congruent set of zs and exponent vectors,
                              ;; x is the product of the zs and y^2 is the sum
                              ;; of the exponent vectors. Get the sqrt of y^2 by halving all
                              ;; the values in the exponent vector.
                              (let [x (->> c-set
                                           (map #(get zs %))
                                           (reduce (fn [a b] (mod-mul a b n))))
                                    y (->> c-set
                                           (mapv #(get exponent-matrix %))
                                           (apply map +)
                                           (map #(quot % 2))
                                           (map (fn [b e] (expt b e)) factor-base)
                                           (reduce (fn [a b] (mod-mul a b n))))]
                                [x y]))]
        (->> (map x-y-pair congruent-sets)
             (mapcat (fn [[x y]] (find-factors x y n)))
             (into (sorted-set)))))))

; 8937589375383974823423423523590825823409839092390489038
; 120778234802486146262478696264740889505538366113384987N
