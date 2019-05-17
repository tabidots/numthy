(ns numthy.modular-arithmetic.modular-fibonacci
  (:require [flatland.useful.seq :refer [take-until]]
            [numthy.primes.is-prime :refer [quick-prime?]]
            [numthy.arithmetic-fns :refer [divisors]]
            [numthy.modular-arithmetic.primitive-roots :refer [powers-of-a-mod-n primitive-root?]]
            [numthy.modular-arithmetic.quadratic-residue :refer [mod-quadratic-zeroes]]
            [numthy.modular-arithmetic.utils :refer [mod-pow]]))

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

;; References
;; https://en.wikipedia.org/wiki/Pisano_period
;; http://webspace.ship.edu/msrenault/fibonacci/fib.htm
;; http://webspace.ship.edu/msrenault/fibonacci/fibfactory.htm
;; https://www.mathstat.dal.ca/FQ/Scanned/15-4/deleon.pdf

(defn modular-fibonacci
  "Lazily generates the first c numbers of the Fibonacci sequence mod n."
  [c n]
  (letfn [(step [[f0 f1]] [(mod f1 n) (mod (+ f0 f1) n)])]
    (take c (map first (iterate step [0 1])))))

(comment
 "A Pisano period is the length of one cycle of the sequence of Fibonacci numbers (mod n).

 e.g. Fibonacci:   [0 1 1 2 3 5 8 13 21 34 55 89 144 233 377 610 987 1597]
      Fib (mod 3): [0 1 1 2 0 2 2 1][0  1  1  2  0   2   2   1   0  ][1..]
      Period:      |-------8-------|

 Interestingly, this sequence is periodic for all n > 1. Though there is no formula to derive
 the period from n, we can see that because every repetition begins with [0 1 ⋯], the
 Pisano period is the smallest integer k > 1 s.t. F_k mod n = 0 and F_{k+1} mod n = 1.")
 ;; Idea adapted from https://medium.com/competitive/huge-fibonacci-number-modulo-m-6b4926a5c836)

(defn pisano-cycle
  "Lazily generates one cycle of the Fibonacci sequence mod n."
  [n]
  (when (> n 1)
    (letfn [(step [[f0 f1]] [(mod f1 n) (mod (+ f0 f1) n)])]
      (cons 0 (map first (take-while #(not= % [0 1]) (iterate step [1 1])))))))

(defn pisano-period
  "The period with which the sequence of Fibonacci numbers (mod n) repeats."
  ;; Using loop-recur is 2-3x faster than lazy sequences here.
  [n]
  (when (> n 1)
    (loop [f0 1 f1 1 k 0] ;; start from F_1 to eliminate checking that k > 1 on every iteration
      (if (and (zero? f0) (= 1 f1)) (inc k) ;; we started 1 index ahead, so we have to inc the result
        (recur (mod f1 n) (mod (+ f0 f1) n) (inc k))))))

(defn- naive-pisano-order
  "The number of zeroes in one Pisano period of n."
  [n]
  (count (filter zero? (pisano-cycle n))))

(defn- naive-pisano-rank
  "The rank of the Fibonacci sequence mod n is the index of its first 0. That is, the least
  positive integer k s.t. F_k ≡ 0 mod n, or alternatively, s.t. n divides F_k."
  [n]
  (/ (pisano-period n) (naive-pisano-order n)))

(defn pisano
  "Computes the period, order, rank, and multiplier of the Pisano cycle mod n and returns
  the results in a map. The period is the length of one Pisano cycle mod n. The order is
  the number of zeros in the cycle. The rank is the index of its first 0. The multiplier
  is the number in the cycle directly after the first 0."
  [n]
  ;; This is messier than the naive functions but much, much faster.
  (when (> n 1)
    (let [m (volatile! {:order 0})]
      (loop [f0 1 f1 1 k 0] ;; start from F_1 to eliminate checking that k > 1 on every iteration
        (if (zero? f0)
          (do
            (vswap! m update :order inc)
            (when-not (:rank @m) (vswap! m assoc :rank (inc k)))
            (when-not (:multiplier @m) (vswap! m assoc :multiplier f1))
            (if (= 1 f1)
              (assoc @m :period (inc k))
              (recur (long (rem f1 n)) (long (rem (+ f0 f1) n)) (inc k)))) ;; rem is more performant than mod
          (recur (long (rem f1 n)) (long (rem (+ f0 f1) n)) (inc k)))))))

(comment
 "A number x is a Fibonacci primitive root mod p, where p is prime, if not only
 do the powers of x mod p generate all values in (ℤ/pℤ,×) (that is, x is a primitive root
 mod p), but x^n + x^(n+1) ≡ x^(n+2) mod p. In other words, the distinct powers of
 x mod p form a Fibonacci sequence.")

(defn- naive-fibonacci-primitive-root?
  "Naïvely determines if x is a Fibonacci primitive root mod p in a lazy, brute-force
  way by offsetting the sequence of powers of x mod p by 1, twice, and testing if
  F_x mod p + F_x+1 mod p = F_x+2 mod p for 0 ≤ x < p."
  [x p]
  (when (primitive-root? x p)
    (let [ps (powers-of-a-mod-n x p)]
      (->> (map #(= (-> (+ %1 %2) (mod p)) %3)
                ps (rest ps) (rest (rest ps)))
           (every? true?)))))

(defn fibonacci-primitive-root?
  "Efficiently determines if x is a Fibonacci primitive root mod p by creating one cycle
  of a Fibonacci sequence starting with [x mod p, x^2 mod p] until it reaches 1, then
  testing that against the sequences of powers of x mod p."
  [x p]
  (letfn [(step [[f0 f1]] [(mod f1 p) (mod (+ f0 f1) p)])]
    (->> (iterate step [x (mod-pow x 2 p)])
         (map first)
         (take-until #(= % 1))
         (= (rest (powers-of-a-mod-n x p))))))

(defn fibonacci-primitive-roots
  "The Fibonacci primitive roots mod p (p prime) are its primitive roots x s.t.
  x^2 - x - 1 ≡ 0 mod p. This is a reverse engineering of the Fibonacci sequence,
  where 1 + x^1 ≡ x^2 mod p. Returns nil if there are no roots."
  [p]
  (when-let [roots (mod-quadratic-zeroes {2 1, 1 -1, 0 -1} p)]
    (not-empty (filter #(primitive-root? % p) roots))))

(defn- naive-has-fib-prim-root?
  "A prime p > 5 has at least one Fibonacci primitive root if p ≡ 1 or 9 (mod 10)
  and its Pisano period is (p - 1). This shortcut minimizes the number of loops to find
  the shortest periodic sequence. Returns true if there are roots, and false/nil otherwise.
  For larger p, this is faster than finding the roots but slower for small p."
  [p]
  (or (= p 5) ;; exception
      (letfn [(kth-fib-pair [k] ;; the kth and k+1th Fib numbers mod p
                            (-> (iterate (fn [[f0 f1]]
                                           (map #(mod % p) [f1 (+ f0 f1)])) [0N 1N])
                                (nth k)))
              (begins-pisano-cycle? [k]
                                    (= [0N 1N] (kth-fib-pair k)))]
        (when (and (quick-prime? p)
                   (contains? #{1 9} (mod p 10))
                   (begins-pisano-cycle? (dec p)))
          ;; Confirm (p-1) is shortest period by testing all divisors d < p of (p-1)
          (not-any? begins-pisano-cycle? (divisors (dec p)))))))

(defn has-fibonacci-primitive-root?
  "Efficiently determines whether p has at least one Fibonacci primitive root.
  Returns nil if p is not prime."
  [p]
  (or (= p 5)
      (and (quick-prime? p)
           (contains? #{1 9} (mod p 10))
           (= (dec p) (pisano-period p)))))
