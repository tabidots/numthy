(ns numthy.modular-arithmetic.utils
  (:require [clojure.math.numeric-tower :refer [gcd expt]]
            [numthy.helpers :refer [pairwise-coprime?]]
            [numthy.primes.is-prime :refer [quick-prime?]]
            [numthy.perfect-powers :refer [perfect-powers]]))

(defn mod-mul
  [a b m]
  (-> (*' a b) (mod m)))

(defn mod-pow
  [b e m]
  (-> (expt b e) (mod m)))

(defn- euclidean-gcd
  "Uses the Euclidean algorithm to find the greatest common divisor of a and b."
  [a b]
  (if (zero? b) a
    (recur b (mod a b))))

(defn- extended-euclidean
  "Given integers a and b, uses the extended Euclidean algorithm to find
  gcd(a, b) and Bézout's coefficients x and y such that ax + by = gcd(a, b)."
  [a b]
  (loop [a a b b s0 1 s 0 t0 0 t 1]
    (if (zero? b) {:gcd a :x s0 :y t0}
      (let [q (quot a b)]
        (recur b (mod a b)
               s (- s0 (* q s))
               t (- t0 (* q t)))))))

(defn mod-inverse
  "Given a modulus m and a coprime integer a, finds the multiplicative
  inverse of a (mod m)—i.e., the integer x such that a*x ≡ 1 (mod m)."
  [a m]
  (let [ee (extended-euclidean a m)]
    (when (= 1 (:gcd ee))
      (let [x (:x ee)]
        (if (neg? x) (+ x m) x)))))

(defn chinese-remainder
  "Solves a system of simultaneous congruences x ≡ a1 mod p1 ≡ a2 mod p2 ...
  when input in the form [[a1 p1] [a2 p2]] ..."
  [cs]
  (let [moduli (map last cs)]
    (when (pairwise-coprime? moduli)
      (let [M (reduce * moduli)]
        (reduce (fn [res [a m]]                        ;; ∑ a*b*b' mod M
                  (let [b (/ M m) b' (mod-inverse b m)]
                    (-> (mod-mul a b M) (mod-mul b' M) (+ res) (mod M))))
                0 cs)))))

(defn bezout-coefs
  "Given integers a and b, returns just Bézout's coefficients x and y."
  [a b]
  (map (extended-euclidean a b) [:x :y]))

(defn odd-prime?
  [n]
  (when (> n 2) (quick-prime? n)))

(defn odd-prime-power?
  "Returns [p k] if n is an odd prime power, k > 0. Returns nil otherwise."
  [n]
  (if (odd-prime? n) [n 1]
    (when-let [p-powers (perfect-powers n)]
      (let [[p k] (first p-powers)]
        (when (odd? p) [p k])))))
