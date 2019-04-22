(ns numthy.modular-arithmetic.quadratic-residue
  (:require [clojure.math.numeric-tower :refer [expt]]
            [numthy.modular-arithmetic.utils :refer [mod-pow mod-inverse odd-prime?]]
            [numthy.helpers :refer [coprime?]]
            [numthy.factorization.core :refer [pollard-factorize]]
            [numthy.polynomials.core :refer [degree]]
            [numthy.p-adic :refer [p-adic-order]]
            [numthy.modular-arithmetic.primitive-roots :refer [primitive-roots]]))

(defn legendre-symbol
  "Legendre symbol for odd prime moduli. Returns 1 if a is a quadratic residue mod p."
  [a p]
  (when (odd-prime?)
    (let [pow (-> (- p 1) (/ 2))  ;; a^((p-1)/ 2) mod p
          sym (mod-pow a pow p)]
      (if (> sym 1) (- sym p)     ;; make all non-residues -1
        sym))))

(defn kronecker-symbol
  "Kronecker symbol of an integer (neg, pos, or zero) a mod n. Returns 1 if a is
  a quadratic residue mod n."
  ;; https://en.wikipedia.org/wiki/Kronecker_symbol
  [a n]
  (if-let [ls (legendre-symbol a n)] ls  ;; for prime moduli p > 2
    (reduce-kv (fn [res p e]
                 (let [x (if (> p 2)
                           (legendre-symbol a p)
                           (nth [0 1 0 -1 0 -1 0 1] (mod a 8)))]
                   (* res (expt x e))))
               (if (neg? a) -1 1)
               (pollard-factorize n))))

(defn perfect-power-of-2?
  ;; Maybe this is useful for something?
  [n]
  (reduce (fn [a b]
            (cond
              (> a n) (reduced false)
              (= a n) (reduced true)
              :else   (* a b)))
          (repeat 2)))

(defn msqrt
  "Uses the Tonelli-Shanks algorithm to find the integer x s.t. x^2 = a (mod p).
  Returns nil if modulus is composite or if there is no solution."
  ;; https://eli.thegreenplace.net/2009/03/07/computing-modular-square-roots-in-python
  ;; http://www.math.vt.edu/people/brown/doc/sqrts.pdf
  [a p]
  (when-let [ls (legendre-symbol a p)]
    (when (and (coprime? a p) (pos? ls))
            ;; Find a number n that is a quadratic nonresidue mod p
      (let [n (->> (drop 2 (range))
                   (filter #(neg? (legendre-symbol % p)))
                   first)
            ;; Represent (p - 1) as (s * 2^e)
            e (p-adic-order 2 (dec p))
            s (/ (dec p) (reduce *' (repeat e 2)))]
        (loop [x (mod-pow a (-> s (+ 1) (/ 2)) p) ;; 1st guess at sqrt
               b (mod-pow a s p)                  ;; 1st guess at fudge factor
               g (mod-pow n s p)
               r e]
                ;; Find the smallest m < r s.t. b^(2^m) ≡ 1 mod p
          (let [m (->> (range r)
                       (filter #(= 1 (mod-pow b (expt 2 %) p)))
                       first)]
            (if (zero? m) x
              (recur
                (-> (mod-pow g (expt 2 (- r m 1)) p)
                    (* x) (mod p))                         ;; x * g^(2^(r-m-1)) % p
                (-> (mod-pow g (expt 2 (- r m)) p)
                    (* b) (mod p))                         ;; b * g^(2^(r-m)) % p
                (-> (mod-pow g (expt 2 (- r m)) p)
                    (mod p))                               ;; g^(2^(r-m)) % p
                m))))))))

(defn mod-quadratic-zeroes
  "Returns the two solutions to any quadratic equation mod p, or nil if none."
  [pnml p]
  (when (and (= 2 (degree pnml)) (.isProbablePrime (biginteger p) 5))
    (let [{a 2, b 1, c 0} pnml
          discriminant    (-> (*' b b)
                              (- (* 4 a c)))]
      (when-let [msr-d (msqrt discriminant p)]
        ;; (-b ± sqrt(4ac)) / 2a ==>> (-b ± msqrt(4ac)) * 2a^-1 mod p
        [(-> (- b) (+ msr-d) (* (mod-inverse (* 2 a) p)) (mod p))
         (-> (- b) (- msr-d) (* (mod-inverse (* 2 a) p)) (mod p))]))))

(defn quadratic-residues
  "Returns all quadratic residues mod m; that is, all integers x < m s.t. there is
  an integer q^2 ≡ x mod m."
  [m]
  (if (odd-prime? m)
    (->> (primitive-roots m)      ;; https://math.stackexchange.com/a/588778
         (map #(mod-pow % 2 m))   ;; if a is p.r. mod p, then a^k is q.r. mod p
         (cons 1)                 ;; where k even.
         (into (sorted-set)))
    (->> (range 1 m)
         (filter #(= 1 (kronecker-symbol % 7)))
         (into (sorted-set)))))
