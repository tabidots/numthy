(ns numthy.factorization.cfrac
  "Implementation of the continued-fraction factorization algorithm, which belongs to the family of congruence-of-squares methods."
  (:require [clojure.core.matrix :as m]
            [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :refer [abs]]
            [numthy.factorization.squares-utils :refer [factorize-rhs find-factor smoothness-bound smooth-factorization]]
            [numthy.factorization.trial-division :refer [batch-smooth-filter]]
            [numthy.helpers :refer [isqrt]]
            [numthy.modular-arithmetic.utils :refer [mod-pow]]
            [numthy.modular-arithmetic.quadratic-residue :refer [quadratic-residue?]]
            [numthy.primes.sieves :refer [primes-to squarefrees-to]]))

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

;; https://yourmaths.wordpress.com/2015/10/05/factoring-large-numbers-with-the-continued-fraction-method/
;; https://programmingpraxis.com/2010/09/14/the-factorization-of-f7-part-1/
;; https://frenchfries.net/paul/factoring/theory/cont.frac.html

(defn- ->small-residue
  "If x^2 mod n ≦ bound, returns x^2 mod n. If n - x^2 mod n ≦ bound, returns -(x^2) mod n.
  Else, returns nil. For the CFRAC algorithm, this bound is 2√(n). "
  [x n bound]
  (let [x2 (mod-pow x 2 n)]
    (if (<= x2 bound) x2
      (let [neg-mod (- x2 n)]
        (when (<= (- bound) neg-mod 0) neg-mod)))))

(defn- find-multiplier
  "Find the k < B that maximizes the number of primes < B of which k*n is a quadratic residue."
  [n B]
  (let [bound (int B)]
    (apply max-key
      (fn [k]
        (->> (primes-to bound)
             (filter (fn [p] (quadratic-residue? (* k n) p)))
             (count)))
      (squarefrees-to bound))))

(defn- sqrt-cfrac-numerators
  "Produces a lazy infinite sequence of numerators of the convergents of
  the continued fraction representation of √n. Note that even though the continued
  fraction is periodic, the numerators are not periodic. They may not be distinct,
  however, when reduced mod n."
  [n]
  (let [a0 (isqrt n)
        cfrac-next (fn [[m d a h0 h1]]
                     (let [m' (-> (*' d a) (-' m))
                           d' (/ (-' n (*' m' m')) d)
                           a' (quot (+' a0 m') d')
                           h' (+' (*' a' h1) h0)]
                       [m' d' a' (mod h1 n) (mod h' n)]))]
    (lazy-seq (map peek (iterate cfrac-next [0 1 a0 1 a0])))))

(defn- add-relations
  "Given a factor base (vector of primes) with B elements, finds at least B+1 values of x from
  the numerators of the convergents of √(n), s.t. x^2 mod n is some value r or -r mod n
  that is [a] ≦ 2√(n) and [b] smooth over the factor base. Those values of x are then factorized
  over the primes in the factor base and returned in a map of the form {x mod n, ̂r mod n}
  where ̂r is the residue of x^2 mod n with the least absolute value."
  [n factor-base limit]
  ;; Originally written using a loop-recur that checked smooths in real-time. It worked
  ;; fast for small numbers but the code was very messy and it was not possible to take advantage of
  ;; the batch smoothness checking algorithm. This code is much cleaner, at the expense of
  ;; a little extra running time on smaller numbers.
  (let [residue-bound (* 2 (isqrt n))]
    (->> (partition 250 (sqrt-cfrac-numerators n))
         (r/reduce (fn [smooths convergents]
                     (if (> (count smooths) limit) (reduced smooths)
                       (->> convergents
                            (batch-smooth-filter factor-base #(->small-residue % n residue-bound))
                            (factorize-rhs factor-base)
                            (into smooths))))
                   {}))))

(defn cfrac
  "Uses the Brillhart-Morrison continued fraction factorization algorithm to find one
  prime factor of n. If the original value of n does not yield a factor, tries again using
  a multiplier k s.t. k*n factors over a maximal number of small primes in the factor base."
  [n]
  (let [B           (smoothness-bound n)
        factor-base (->> (primes-to (inc (int B)))
                         (filter #(quadratic-residue? n %))
                         (cons -1))
        quota       (inc (count factor-base))]
    (or (->> (add-relations n factor-base quota)
             (find-factor n factor-base))
        (let [k (find-multiplier n B)]
          (->> (add-relations (* k n) factor-base quota)
               (find-factor n factor-base))))))
