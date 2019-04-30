(ns numthy.factorization.cfrac
  "Implementation of the continued-fraction factorization algorithm, which belongs to the family of congruence-of-squares methods."
  (:require [clojure.math.numeric-tower :refer [abs expt]]
            [numthy.factorization.squares-utils :refer [smooth? smoothness-bound find-factor-from-congruent-relations]]
            [numthy.helpers :refer [isqrt]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.modular-arithmetic.utils :refer [mod-pow]]
            [numthy.modular-arithmetic.quadratic-residue :refer [quadratic-residue?]]))

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

;; https://yourmaths.wordpress.com/2015/10/05/factoring-large-numbers-with-the-continued-fraction-method/
;; https://programmingpraxis.com/2010/09/14/the-factorization-of-f7-part-1/
;; https://frenchfries.net/paul/factoring/theory/cont.frac.html

(defn- ->small-residue
  "Returns x^2 if |x^2| mod n is within some bound, and nil otherwise.
  For the CFRAC algorithm, this bound is 2√(n). "
  [x n bound]
  ;; According to this paper (https://www.rose-hulman.edu/~bryan/lottamath/factor.pdf)
  ;; negative residues mod n should be included, and -1 should also be included in the factor
  ;; base. However, this resulted in slightly worse performance, not better.
  (let [x2 (mod-pow x 2 n)]
    (if (<= x2 bound) x2
      (let [neg-mod (abs (- x2 n))]
        (when (<= neg-mod bound) neg-mod)))))

(defn- generate-relations
  "Given a factor base (vector of primes) with B elements, finds B+1 values of q from
  the numerators of the convergents of √(n), s.t. q^2 mod n is some value r or -r mod n
  that is [a] ≦ 2√(n) and [b] smooth over the factor base."
  [n factor-base num-relations]
  ;; This function is pretty ugly because it needs to generate convegents and check their
  ;; eligibility in parallel. Doing this sequentially would be much cleaner, but far less efficient
  ;; because we don't know how many convergents it will take in order to yield a given number
  ;; of eligible ones.
  ;; Algo reference: https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html
  (let [prime-map (into (sorted-map) (zipmap factor-base (repeat 0)))
        x         (isqrt n)
        two-x     (* 2 x)]
    (loop [y x
           z 1N
           a two-x
           e1 1N e2 0N
           f1 0N f2 1N
           relations (transient {})]
      (let [yy (-> (*' a z) (-' y))
            zz (-> (-' n (*' yy yy)) (/ z) biginteger)
            aa (-> (+' x yy) (/ zz) biginteger)
            h  (-> (+' e2 (*' x f2)) (mod n))]   ;; numerator of convergent
        (if (or (= zz 1) (> (count relations) num-relations))
          (do (println (count relations) "relations collected.")
            (persistent! relations))
          (recur
            yy zz aa
            e2 (-> (*' aa e2) (+' e1) (mod n))
            f2 (-> (*' aa f2) (+' f1) (mod n))
            (if-let [expv (some-> (->small-residue h n two-x)
                                  (smooth? factor-base prime-map))]
              (assoc! relations h expv)
              relations)))))))

(defn- cfrac-single
  "Uses the CFRAC factorization algorithm to find at least one factor of n.
  Returns nil if algorithm fails to find a factor with the given multiplier.
  Default multiplier is 1."
  [N multiplier]
  (let [kN          (* N multiplier)
        B           (smoothness-bound kN)
        factor-base (filter #(and (prime? %) (quadratic-residue? kN %)) (range B))
        ;pool-size   (expt (Math/log10 n) 2)  This is a fail-safe upper bound but too large for large n and thus too slow
        ;_           (println "Looking for" (inc (long pool-size)) "relations...")
        relations   (generate-relations kN factor-base (inc (count factor-base)))]
    (find-factor-from-congruent-relations N factor-base relations)))

(defn cfrac-factorize
  [N]
  (some #(cfrac-single N %) (range 1 20)))

;; TODO: How to definitively find the best multiplier k?
