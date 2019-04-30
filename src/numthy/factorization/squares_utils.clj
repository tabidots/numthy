(ns numthy.factorization.squares-utils
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
  [n prime-map]
  (loop [x      n
         m      prime-map
         primes (keys m)]
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
  ;; ↓ The idxs of leading 1s tell us the idxs of linearly independent
  ;; rows in the exponent matrix.
  [row]
  (first (keep-indexed (fn [i x]
                         (when (= 1 x) i))
                       row)))

(defn linearly-independent-idxs
  [m]
  (->> (m/transpose m)
       (linalg/reduced-row-echelon)
       (r/map index-of-leading-1)
       (r/remove nil?)
       (into [])))

(defn linearly-independent-idxs-gf2
  [m]
  (->> (m/transpose m)
       (linalg/rref-gf2)
       (r/map index-of-leading-1)
       (r/remove nil?)
       (into [])))

(defn make-x-from-idxs
  [idxs zs n]
  (->> (m/order zs idxs)
       (r/fold (fn ([] 1) ([a b] (mod-mul a b n))))))

(defn make-y-from-idxs
  [idxs expm factor-base n]
  (->> (m/order expm idxs)
       (apply m/add)
       (zipmap factor-base)
       (r/fold (r/monoid (fn [a b] (mod-mul a b n)) (constantly 1))
               (fn ([res b e]
                    (-> (expt b (quot e 2))  ;; Adding the exponent vectors gives y^2
                        (mod-mul res n)))))))

(defn find-factor-from-congruent-relations
  "Segregates the rows of an exponent matrix, reduced mod 2, into:
  - A hash-map of its linearly independent rows (in the form {idx [exp-vec]})
  - A hash-map of its linearly dependent rows (in the form {idx [exp-vec]})
  This redundancy is due to the fact that the original, non mod-reduced matrix
  is also needed to find the congruences. Given a matrix segregated into linearly independent and dependent rows,
  where A is the matrix of linearly independent rows, and b is a dependent row,
  this solves Ax = b for each b, mod 2. x is a vector containing 0s and 1s,
  where 1s mean keep the corresponding row in A, and 0s mean to drop it.
  Returns the idxs of the filtered rows in A with idx of b appended. generate pairs of x and y such that x^2 ≡ y^2 mod n,
   given a congruent set of zs and exponent vectors. x is the product of the zs
   and y^2 is the sum of the exponent vectors. Get the sqrt of y^2 by halving all
   the values in the exponent vector. Given x, y, and modulus n, uses congruence of squares mod n to determine whether
  gcd(x + y, n) and gcd(x - y, n) are factors of n. This function does not exclude
  trivial factors (1 and n)"
  [n factor-base relations]
  (let [zs               (vec (keys relations))
        expm             (vec (vals relations))
        binary-expm      (mapv row->mod2 expm)
        ;; ↓ "To find the linearly independent rows of the exponent matrix,"
        ;; "we can get the RREF of the exponent matrix transposed"
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
              ;; Winners are a solution to the matrix equation Ax = b and look like [0 1 1 0 1]
              idxs (->> (linalg/solve-gf2 At b) ; old: (row->mod2 (linalg/least-squares At b))
                        ;; If the independent indxs are [1 2 3 4 5] and winners are [0 1 1 0 1],
                        ;; then congruent indxs are [- 2 3 - 5] -> [2 3 5].
                        (keep-indexed (fn [i x] (when-not (zero? x)
                                                  (get independent-idxs i))))
                        ; Add the index for row b. e.g. (10 2 3 5) if solving for row 10
                        (cons d-idx)
                        (vec))
              x    (make-x-from-idxs idxs zs n)
              y    (make-y-from-idxs idxs expm factor-base n)
              fact (gcd (+ x y) n)]
          (if (and (not= (mod x n) (mod y n)) (< 1 fact n)) fact
            (recur (rest dependent-idxs))))
        "No factors found"))))
