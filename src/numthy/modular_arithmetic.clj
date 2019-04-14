(ns numthy.modular-arithmetic
  (:require [clojure.math.numeric-tower :as tower]
            [clojure.math.combinatorics :as combo]
            [numthy.helpers :as h]))

;; BASIC OPERATIONS
;; Modular multiplication and exponentiation basically use the same algorithm,
;; with the only differences being the operation (+ or *) and initial result (0 or 1).

(defn mod-mul
  ;; Translated from https://www.geeksforgeeks.org/multiply-large-integers-under-large-modulo/
  "Quickly calculates a * b % m. Useful when a and b are very large integers."
  [a b m]
  (if (or (= m 1) (= a m)) 0
    (loop [a (mod a m) b b res 0]
      (if (zero? b) res
        (recur (-> (+' a a) (mod m))
               (.shiftRight (biginteger b) 1)
               (if (odd? b)
                 (-> (+' res a) (mod m))
                 res))))))

(defn mod-pow
  ;; Adapted from https://en.wikipedia.org/wiki/Modular_exponentiation
  "Quickly calculates a ^ b % m. Useful when a and b are very large integers."
  [base exp m]
  (if (or (= m 1) (= base m)) 0
    (loop [base (mod base m) e exp res 1]
      (if (zero? e) res
        (recur (-> (*' base base) (mod m))
               (.shiftRight (biginteger e) 1)
               (if (odd? e)
                 (-> (*' res base) (mod m))
                 res))))))

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

;; GROUPS AND GENERATORS

(defn complete-residue-system?
  "A complete residue system mod m is a set containing precisely one
  representative of each residue class (congruence class) mod m.
  i.e., {-13 4 17 18} mod 4 = {3 0 1 2} ∴ true"
  [coll m]
  (= (set (range m)) (set (map #(mod % m) coll))))

(defn generators-additive ;; https://www.doc.ic.ac.uk/~mrh/330tutor/ch06.html
  "Generators of (ℤ/nℤ,+), the additive group of integers mod n—
  i.e., the integers less than n that generate all integers < n
  when added repeatedly to themselves, excluding 0 (universal non-generator) and
  1 (universal generator)."
  [n]
  (let [group    (set (range n))
        generate (fn [seed] (set (map #(mod-mul seed % n) (range n))))]
    (->> (drop 2 (range n))
         (filter #(= group (generate %))))))

(defn coprime? [a b]
  (= 1 (tower/gcd a b)))

(defn multiplicative-group
  "(ℤ/nℤ,×), the multiplicative group of integers mod n,
  or all the integers < n that have a multiplicative inverse mod n.
  Effectively, it is all of the integers less than n that are coprime to n."
  [n]
  (if (h/prime? n)
    (range 1 n)
    (filter (partial coprime? n) (range 1 n))))

(defn powers-of-a-mod-n
  "a^k (mod n) for all 0 ≦ k < n, where k ∈ ℤ."
  [a n]
  (map #(mod-pow a % n) (range n)))

(defn naive-primitive-roots
  "Primitive roots of (ℤ/nℤ,×), the multiplicative group of integers mod n,
  i.e., all of the integers < n that, when multiplied repeatedly by themselves,
  generate all integers < n. This is essentially a brute-force method
  that will take a long time for n > 500 or so."
  [n]
  (let [group (set (multiplicative-group n))
        roots (filter #(= group (set (powers-of-a-mod-n % n))) group)]
    (when (not-empty roots)
      (into (sorted-set) roots))))

(defn primitive-root?
  "Efficiently determines if x is a primitive root mod p, where p is prime,
  by leveraging the fact that x is a primitive root mod p if no prime factor of
  Φ(p) (i.e., p - 1) is a primitive root mod p."
  ; https://stackoverflow.com/a/26636457/4210855
  [x p]
  (let [phi (dec p)
        pfs (filter h/prime? (h/factors phi))]
    (not-any? #(= 1 (mod-pow x (/ phi %) p)) pfs)))

(defn cyclic?
  "Returns true if (ℤ/nℤ,×) ≅ C_Φ(n), i.e., the multiplicative group of
  integers mod n is cyclic. That is, the group has at least one primitive root.
  Returns nil otherwise."
  [n]
  (when (.isProbablePrime (biginteger n) 5)
    (some #(primitive-root? % n) (range 1 n))))

(defn multiplicative-order
  "ord_n(a), the smallest positive integer k such that a^k ≡ 1 (mod n) where
  a is coprime to n. a^0 ≡ 1 (mod n) for any n. "
  [a n]
  (when (coprime? a n)
    (first (filter #(= 1 (mod-pow a % n)) (rest (range))))))

(comment
  "If we know we are dealing with ℤ/pℤ, this is faster."
  (first (filter #(= 1 (mod-pow a % p)) (pollard-factorize (dec p)))))

(comment
  "For multiplicative order, a quick way to find k for small n is to count
  the distinct powers of a (mod n). However, this is intractable for large n."
  (count (distinct (powers-of-a-mod-n a n))))

(comment
  "Multiplicative order also enables an alternative test for primitive root-ness.
  Works for prime and composite moduli, but will be impractically slow for large moduli."
  (= (multiplicative-order x m) (h/phi p)))

(defn pairwise-coprime?
  "True if every distinct pair of (a, b) in xs is coprime."
  [xs]
  (every? (partial apply coprime?) (combo/combinations xs 2)))

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

(defn pollard-rho
  "Uses Pollard's rho algorithm to factorize a large integer n."
  [n]
  (letfn [(g [x] (+ 3 (mod-pow x 2 n)))]    ;; g(x) = (x^2 + 3) mod n
    (loop [a (g 2) b (g (g 2)) d 1]
      (cond
        (= d n)   nil
        (< 1 d n) d
        :else     (recur
                    (g a) (g (g b)) (tower/gcd (- a b) n))))))

(defn pollard-factorize
  "Uses Pollard's rho algorithm recursively to find many (but not necessarily all)
  factors of a large integer n. Much faster than trial division for n > 10 million."
  ;; Try it with 45883247897, or 1240087781 (a prime)
  [n]
  (loop [n   n
         res (conj (sorted-set) 1 n)] ;
    (if (< n 10000000)
      (into res (h/factors n))
      (if-let [p (pollard-rho n)]
        (let [q (/ n p)] ;; n is composite
          (recur (max p q) (conj res p q)))
        (conj res n))))) ;; n is prime

(defn carmichael
  "The smallest positive integer k such that a^k ≡ 1 (mod n) for all a coprime to n."
  [n]
  (if (= n 1) 1
    (letfn [(prime-power-carm [[prime power]]
              ;; for a power of an odd prime and for 2 and 4, λ(n) is equal to
              ;; the Euler totient φ(n); for powers of 2 greater than 4 it is
              ;; equal to half of the Euler totient.
              (cond
                (and (odd? prime)
                     (= 1 power)) (dec prime) ;; speed up calculation
                (odd? prime)      (h/phi (tower/expt prime power))
                (> power 2)       (/ (h/phi (tower/expt prime power)) 2)
                :else             (h/phi (tower/expt prime power))))]
      ;; λ(n) is the least common multiple of the λ of each of its prime power factors
      (->> (h/prime-factorization n)
           frequencies
           (map prime-power-carm)
           (reduce tower/lcm)))))

;; TODO: discrete logarithm problem
;; https://www.doc.ic.ac.uk/~mrh/330tutor/ch06s02.html

;; TODO: Brent's cycle finding algorithm
;; TODO: Pollard-Brent factorization https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf

;; TODO: https://en.wikipedia.org/wiki/Carmichael_number
;; TODO: https://en.wikipedia.org/wiki/Euler_pseudoprime
;; https://en.wikipedia.org/wiki/Fermat%27s_little_theorem

;; https://en.wikipedia.org/wiki/Congruent_number

;; TODO: Dirichlet character https://en.wikipedia.org/wiki/Dirichlet_character

(defn bezout-coefs
  "Given integers a and b, returns just Bézout's coefficients x and y."
  [a b]
  (map (extended-euclidean a b) [:x :y]))
