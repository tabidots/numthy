(ns numthy.factorization.cfrac
  "Implementation of the continued-fraction factorization algorithm, which belongs to the family of congruence-of-squares methods."
  (:require [clojure.math.numeric-tower :refer [abs expt gcd]]
            [clojure.core.matrix :as m]
            [clojure.core.reducers :as r]
            [numthy.factorization.squares-utils :refer [smooth? smoothness-bound find-factor]]
            [numthy.helpers :refer [isqrt divisible?]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.primes.sieves :refer [primes-to squarefrees-to]]
            [numthy.modular-arithmetic.utils :refer [mod-pow mod-mul]]
            [numthy.modular-arithmetic.quadratic-residue :refer [quadratic-residue? kronecker-symbol]]
            [numthy.factorization.trial-division :refer [batch-smooth-filter]]))

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
            zz (-> (-' n (*' yy yy)) (/ z) bigint)
            aa (-> (+' x yy) (/ zz) bigint)
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

(defn find-multiplier
  "Find the k < B that maximizes the number of primes < B of which k*n is a quadratic residue."
  [n B]
  (let [bound (int B)]
    (apply max-key
      (fn [k]
        (->> (primes-to bound)
             (filter (fn [p] (= 1 (kronecker-symbol (* k n) p))))
             (count)))
      (squarefrees-to bound))))

(defn cfrac
  [n]
  (let [B           (smoothness-bound n)
        factor-base (cons -1 (filter #(and (prime? %) (quadratic-residue? n %)) (range B)))
        _           (println "Finding relations (k = 1)...")
        relations   (generate-relations n factor-base (inc (count factor-base)))]
    (if-some [factor (find-factor n factor-base relations)]
      factor
      (let [k (find-multiplier n B)
            _ (println "Finding relations (k =" k ")...")
            relations (generate-relations (* k n) factor-base (inc (count factor-base)))]
        (find-factor n factor-base relations)))))

(comment
 "Below is experimental code, not fully working yet. Batch smoothness checking still not working correctly"
 (defn smooth-factorization
   [n primes]
   (let [n' (abs n)]
     (apply merge (for [p primes]
                    (loop [x n' exp 0]
                      (cond
                        (neg? p)         (if (neg? n) {-1 1} {-1 0})
                        (pos? (mod x p)) {p exp}
                        :else            (recur (/ x p) (inc exp))))))))

 (defn- batch-convergents
   "Given a factor base (vector of primes) with B elements, finds B+1 values of q from
  the numerators of the convergents of √(n), s.t. q^2 mod n is some value r or -r mod n
  that is [a] ≦ 2√(n) and [b] smooth over the factor base."
   [n batch-size]
   ;; This function is pretty ugly because it needs to generate convegents and check their
   ;; eligibility in parallel. Doing this sequentially would be much cleaner, but far less efficient
   ;; because we don't know how many convergents it will take in order to yield a given number
   ;; of eligible ones.
   ;; Algo reference: https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html
   (let [x         (isqrt n)
         two-x     (* 2 x)]
     (loop [y x
            z 1N
            a two-x
            e1 1N e2 0N
            f1 0N f2 1N
            hs (transient [])]
       (let [yy (-> (*' a z) (-' y))
             zz (-> (-' n (*' yy yy)) (/ z) bigint)
             aa (-> (+' x yy) (/ zz) bigint)]
         (if (>= (count hs) batch-size) (persistent! hs)
           (recur
             yy zz aa
             e2 (-> (*' aa e2) (+' e1) (mod n))
             f2 (-> (*' aa f2) (+' f1) (mod n))
             (conj! hs (-> (+' e2 (*' x f2)) (mod n)))))))))

 (defn batch-relations
   [n factor-base]
   (let [batch-size    (expt (Math/log10 n) 2.9)
         residue-bound (*' 2 (isqrt n))
         hs            (batch-convergents n batch-size)
         _ (println (count hs) "convergents collected.")
         smooths       (batch-smooth-filter factor-base
                                            #(->small-residue-signed % n residue-bound)
                                            hs)
         _ (println (count smooths) "smooths found.")]
      (when (not-empty smooths)
        (reduce-kv (fn [relations x2 x2modn]
                     (assoc relations x2 (->> (smooth-factorization x2modn factor-base)
                                              vals
                                              vec)))
                   {} smooths))))

 (time (let [n 16850989]
         (let [B             (smoothness-bound n)
               factor-base   (cons -1 (filter #(and (prime? %) (quadratic-residue? n %)) (range B)))]
           (count (batch-relations n factor-base))))))
