(ns numthy.primes.generate)

(defn rand-prime
  "Use Java interop to generate a random prime of the specified number of bits."
  [bits]
  (BigInteger/probablePrime bits (java.util.Random.)))

(defn rand-safe-prime
  "Use Java interop to generate a random prime of the specified number of bits,
  that is also of the form 2p + 1, where p is also prime."
  [bits]
  (->> #(rand-prime bits)
       repeatedly
       (filter #(.isProbablePrime (biginteger (/ (dec %) 2)) 5))
       first))
