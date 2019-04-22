(ns numthy.primes.generate)

(defn rand-prime
  "Use Java interop to generate a random prime of the specified number of bits."
  [bits]
  (BigInteger/probablePrime bits (java.util.Random.)))
