(ns numthy.modular-arithmetic.modular-fibonacci
  (:require [numthy.primes.is-prime :refer [prime?]]
            [numthy.arithmetic-fns :refer [divisors]]
            [numthy.modular-arithmetic.primitive-roots :refer [powers-of-a-mod-n primitive-root?]]
            [numthy.modular-arithmetic.quadratic-residue :as qr]))

;; First, some naïve implementations to explore the concepts in a more
;; understandable way. These will only work in reasonable time for small inputs.
;; More practical implementations follow.

(defn fibonacci-primitive-root?
  "Naïvely determines if x is a Fibonacci primitive root mod p, where p is prime,
  which means that not only do the powers of x mod p generate all values in
  (ℤ/pℤ,×), but x^n + x^(n+1) ≡ x^(n+2) mod p."
  [x p]
  (when (primitive-root? x p)
    (let [ps (powers-of-a-mod-n x p)]
      (->> (map #(= (-> (+ %1 %2) (mod p)) %3)
                ps (rest ps) (rest (rest ps)))
           (every? true?)))))

(defn pisano-period
  "The period with which the sequence of Fibonacci numbers (mod n) repeats."
  ;; e.g. Fibonacci:   [0 1 1 2 3 5 8 13 21 34 55 89 144 233 377 610 987 1597]
  ;;      Fib (mod 3): [0 1 1 2 0 2 2 1][0  1  1  2  0   2   2   1   0  ][1..]
  ;;      Period:      |-------8-------|
  ;; This sequence is periodic for all n, and because every repetition begins with
  ;; [0 1 ⋯], the Pisano period is also the smallest integer k > 1 s.t. the
  ;; k-th Fib-num mod n is 0 and the k+1-th Fib-num mod n is 1.
  ;; Shortcut idea taken from https://medium.com/competitive/huge-fibonacci-number-modulo-m-6b4926a5c836
  [n]
  (when (> n 1)
    (letfn [(next-pair [[f0 f1]] [f1 (+ f0 f1)])]
      (loop [fp [0N 1N] k 0]
        (if (and (> k 1) (= (map #(mod % n) fp) [0N 1N])) k
          (recur (next-pair fp) (inc k)))))))

(defn naive-has-fib-prim-root?
  "A prime p > 5 has at least one Fibonacci primitive root if it is ≡ 1 or 9 (mod 10)
  and its Pisano period is (p - 1). This shortcut minimizes the number of loops to find
  the shortest periodic sequence."
  [p]
  (if (= p 5) true ;; exception
    (letfn [(kth-fib-pair [k] ;; the kth and k+1th Fib numbers mod p
              (-> (iterate (fn [[f0 f1]]
                             (map #(mod % p) [f1 (+ f0 f1)])) [0N 1N])
                  (nth k)))
            (begins-pisano-cycle? [k]
              (= [0N 1N] (kth-fib-pair k)))]
      (when (and (prime? p)
                 (contains? #{1 9} (mod p 10))
                 (begins-pisano-cycle? (dec p)))
        ;; Confirm (p-1) is shortest period by testing all factors of (p-1)
        (not-any? begins-pisano-cycle? (drop-last (divisors (dec p))))))))

;; TODO: https://en.wikipedia.org/wiki/Pisano_period#Generalizations

(defn has-fib-prim-root?
  [p]
  (when-let [roots (qr/mod-quadratic-zeroes {2 1, 1 -1, 0 -1} p)]
    (some #(primitive-root? % p) roots)))
