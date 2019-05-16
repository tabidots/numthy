(ns numthy.factorization.cfrac
  "Implementation of the continued-fraction factorization algorithm, which belongs to the family of congruence-of-squares methods."
  (:require [clojure.core.matrix :as m]
            [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :refer [abs expt]]
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
  "If x^2 mod n ≤ bound, returns x^2 mod n. If n - x^2 mod n ≤ bound, returns -(x^2) mod n.
  Else, returns nil. For the CFRAC algorithm, this bound is 2√(n). "
  [x n bound]
  (let [x2 (mod-pow x 2 n)]
    (if (<= x2 bound) x2
      (let [neg-mod (- x2 n)]
        (when (<= (- bound) neg-mod 0) neg-mod)))))

(defn- sqrt-cfrac-numerators
  "Produces a lazy sequence of numerators of the convergents of one period of the continued fraction
  representation of √n, reduced mod n."
  [n]
  (let [a0 (isqrt n)
        cfrac-next (fn [[m d a h0 h1]]
                     (let [m' (-> (*' d a) (-' m))
                           d' (/ (-' n (*' m' m')) d)
                           a' (quot (+' a0 m') d')
                           h' (+' (*' a' h1) h0)]
                       [m' d' a' (mod h1 n) (mod h' n)]))]
    (->> [0 1 a0 1 a0]
         (iterate cfrac-next)
         (rest)                         ;; drop the 0th iteration
         (take-while #(> (second %) 1)) ;; truncate to one period of the c-frac (d = 1 marks start of cycle, incl. 0th iter)
         (map peek)
         (lazy-seq))))

(defn- add-relations
  "Given a factor base (vector of primes) with B elements, finds at least B+1 values of x from
  the numerators of the convergents of √(n), s.t. x^2 mod n is some value r or -r mod n
  that is [a] ≤ 2√(n) and [b] smooth over the factor base. Those values of x are then factorized
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
    (some #(->> (add-relations (*' % % n) factor-base quota)
                (find-factor n factor-base))
          (rest (range)))))

(comment
 "The following functions were used to refine the multiplier search. I have left them here for reference.
 The papers all say that the best multipliers are square-free, but I have found that in a minority of
 cases, square-free multipliers do not find a factor, and even in the majority of cases where they do,
 the so-called 'optimal k' that maximizes the number of primes < B of which k*n is a quadratic residue
 does not necessarily find a factor of n the fastest.

 Oddly enough, a factor is always found fastest when k is a perfect square. Go figure!"

 (defn- sqfree-multiplier
   "Find the squarefree integer k < B that maximizes the number of primes < B of which k*n is a quadratic residue."
   [n B]
   (let [bound (int B)]
     (apply max-key
       (fn [k]
         (->> (primes-to bound)
              (filter (fn [p] (quadratic-residue? (*' k n) p)))
              (count)))
       (squarefrees-to bound))))

 (defn- find-multiplier
   "Find the k < log_10(n) that maximizes the number of primes < B of which k*n is a quadratic residue."
   [n B]
   (let [bound (int B)]
     (apply max-key
       (fn [k]
         (->> (primes-to bound)
              (filter (fn [p] (quadratic-residue? (*' k n) p)))
              (count)))
       (range (Math/log10 n)))))

 (defn- test-multipliers
   [n]
   (do
     (println (format "Finding the best k-multiplier for %d..." (biginteger n)))
     (let [B           (smoothness-bound n)
           factor-base (->> (primes-to (inc (int B)))
                            (filter #(quadratic-residue? n %))
                            (cons -1))
           quota       (inc (count factor-base))
           little-k    (find-multiplier n B)
           big-k       (sqfree-multiplier n B)
           squares     (->> (map #(*' % %) (rest (range)))
                            (take-while #(< % big-k))
                            vec)
           mults       (conj squares little-k big-k)]
       (doseq [k mults]
         (println (format "k = %d" (biginteger k)))
         (let [f (volatile! nil)
               s (->> (add-relations (*' k n) factor-base quota)
                      (find-factor n factor-base)
                      (vreset! f)
                      time
                      with-out-str)]
           (if (some? @f)
             (println (format "Factor found! Took %s ms." (re-find #"\d[\S]*" s)))
             (println "No factor found. :("))
           (println "")))))))
