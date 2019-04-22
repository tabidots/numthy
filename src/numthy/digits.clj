(ns numthy.digits)

(defn digit-sum
  "Sum of digits of an integer x in base b. Defaults to base 10."
  ([x]
   (digit-sum x 10))
  ([x b]
   (loop [x x sum 0]
     (if (< x b) (+ sum x)
       (recur (quot x b) (+ sum (rem x b)))))))

(defn digit-product
  "Product of digits of an integer x in base b. Defaults to base 10."
  ([x]
   (digit-product x 10))
  ([x b]
   (loop [x x prod 1]
     (if (< x b) (*' prod x)
       (recur (quot x b) (*' prod (rem x b)))))))

(defn digital-root
  "Quick way, using modular congruences, to compute the digital root of a number,
  which is the single-digit number eventually reached after repeatedly taking its
  digital sum."
  [n]
  (cond
    (zero? n)         0
    (zero? (mod n 9)) 9
    :else             (mod n 9)))

(defn additive-persistence
  "The additive persistence of an integer is the number of times its digital sum
  must be taken to reduce it to a single-digit number."
  [x]
  (loop [x x pers 0]
    (if (< x 10) pers
      (recur (digit-sum x) (inc pers)))))

(defn- loop-mult-d-roots
  [x]
  (loop [x x pers 0]
    (if (< x 10) {:mdr x :persistence pers}
      (recur (digit-product x) (inc pers)))))

(defn multiplicative-persistence
  "The multiplicative persistence of an integer is the number of times its digital
  product must be taken to reduce it to a single-digit number."
  [x]
  (:persistence (loop-mult-d-roots x)))

(defn multiplicative-digital-root
  "The multiplicative digital root of a number is the single-digit number eventually
  reached after repeatedly taking its digital product."
  [x]
  (:mdr (loop-mult-d-roots x)))
