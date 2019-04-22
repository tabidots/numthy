(ns numthy.pseudorandom
  (:require [numthy.primes.generate :refer [rand-prime]]
            [numthy.modular-arithmetic.utils :refer [mod-pow]]
            [numthy.modular-arithmetic.primitive-roots :refer [primitive-root?]]))

;; https://en.wikipedia.org/wiki/Category:Pseudorandom_number_generators
;; https://en.wikipedia.org/wiki/Category:Cryptographically_secure_pseudorandom_number_generators

(defn bin->dec
  "Converts a vector of bits to an integer."
  [bit-vec]
  (reduce (fn [res bit]
            (+' (*' 2 res) bit))
          0 bit-vec))

(defn blum-micali
  "Uses the Blum-Micali algorithm to generate a cryptographically secure pseudorandom number,
  or perhaps more than pseudorandom, since this function generates pseudorandom initial
  values for the algorithm as well."
  [iters]
  (let [p      (rand-prime 16)
        output #(if (<= % (/ (dec p) 2)) 1 0)
        g      (first (filter #(primitive-root? % p) (range)))]
    (loop [x (rand-int 10000) bits []]
      (if (= (count bits) iters) (bin->dec bits)
        (let [a (mod-pow g x p)]
          (recur a (conj bits (output a))))))))
