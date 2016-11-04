(ns hidden-watermark.core
  (:import [java.awt Image]
           [java.awt.image BufferedImage]
           [edu.emory.mathcs.jtransforms.fft DoubleFFT_2D RealFFTUtils_2D])
  (:require [clojure.java.io :as io]
            [mikera.image.core :as img]
            [clojure.core.matrix :as mat])
  (:use  [clojure.core.matrix]
         [clojure.core.matrix.operators]))

(mat/set-current-implementation :vectorz)

(def path "test.png")
(def wm-path "wm.png")

(defn buffer-image [path]
  (let[file (io/file path)]
    (img/load-image file)))

(defn real ^double [^double r ^double i] r)
(defn imaginary ^double [^double r ^double i] i)

(defn dft-fn [ri-fn M]
    (let [w (mat/column-count M)
          h (mat/row-count M)
          ds (mat/to-double-array M)
          ^doubles result (make-array Double/TYPE (* 2 (count ds)))
          dfft-2d (DoubleFFT_2D. h w)
          unpacker (RealFFTUtils_2D. h w)
          nM (mat/new-matrix h w)]
      (System/arraycopy ds 0 result 0 (count ds))
      (.realForward dfft-2d result)
      (dotimes [y h]
        (dotimes [x w]
          (let [r (.unpack unpacker y (* 2 x) result 0)
                i (.unpack unpacker y (inc (* 2 x)) result 0)]
            (mat/mset! nM y x (ri-fn r i)))))
      nM))

(def dft-reals (partial dft-fn real))
(def dft-imaginarys (partial dft-fn imaginary))

;; R is real matrix
;; I is imaginary matrix
(defn idft-ri-fn [ri-fn R I]
  (let [w (mat/column-count R)
        hw (int (/ w 2))
        h (mat/row-count R)
        hh (int (/ h 2))
        ^doubles result (make-array Double/TYPE (* 2 w h))
        dfft-2d (DoubleFFT_2D. h w)
        nM (mat/new-matrix h w)]
    (dotimes [y h]
      (dotimes [x h]
        (aset result (+ (* y 2 h) (* x 2)) ^double (mat/mget R y x))
        (aset result (+ (* y 2 h) (inc (* x 2))) ^double (mat/mget I y x))))
    (.complexInverse dfft-2d result true)
    (dotimes [y h]
      (dotimes [x w]
        (let [r (aget result (+ (* y 2 h) (* x 2)))
              i (aget result (+ (* y 2 h) (inc (* x 2))))]
          (mat/mset! nM y x (ri-fn r i)))))
    nM))

(def idft-ri-reals (partial idft-ri-fn real))

(defn image-to-matrix [^BufferedImage bi rgb-to-val-fn]
  (let [h (.getHeight bi)
        w (.getWidth bi)
        M (mat/new-matrix h w)]
    (dotimes [y h]
      (dotimes [x w]
        (mat/mset! M y x (rgb-to-val-fn (.getRGB bi x y)))))
    M))

(defn matrix-to-image [M val-to-rgb-fn]
  (let [w (mat/column-count M)
        h (mat/row-count M)
        bi (BufferedImage. w h BufferedImage/TYPE_INT_ARGB)]
    (dotimes [y h]
      (dotimes [x w]
        (let [v (mat/mget M y x)]
          (.setRGB bi x y (unchecked-int (val-to-rgb-fn v))))))
    bi))

;; encode here
(def test-matrix (image-to-matrix (buffer-image path) #(-> %)))

;; (def fft-matrix-m (dft-magnitudes test-matrix))
;; (def fft-matrix-p (dft-phases test-matrix))
(def fft-matrix-r (dft-reals test-matrix))
(def fft-matrix-i (dft-imaginarys test-matrix))

(def ifft-matrix (idft-ri-reals fft-matrix-r fft-matrix-i))

;; (img/show (matrix-to-image fft-matrix-r (fn [x] x)) :title "real fft image")
;; (img/show (matrix-to-image fft-matrix-i (fn [x] x)) :title "imaginary fft image")
;; (img/show (matrix-to-image fft-matrix-m (fn [x] x)) :title "magnitude fft image")
;; (img/show (matrix-to-image fft-matrix-p (fn [x] x)) :title "phase fft image")
;; (img/show (matrix-to-image ifft-matrix (fn [x] x)) :title "inverse fft image")

(defn update-submatrix [origin submatrix [x y]]
  (let [origin-dup (mat/clone origin)]
    (mat/assign! (mat/submatrix origin-dup x (mat/row-count submatrix) y (mat/column-count submatrix)) submatrix)
    origin-dup))

(defn display-from-matrix [m title]
  (let [image (matrix-to-image m #(-> %))]
     (img/show image :title title)))

(display-from-matrix fft-matrix-r "real fft")
(display-from-matrix fft-matrix-i "imaginary fft")
(display-from-matrix ifft-matrix "ifft origin image")

(def scaled-wm (img/resize (buffer-image wm-path) 256 256))
(def wm-matrix (image-to-matrix scaled-wm #(-> (/ % 1000000))))

(display-from-matrix wm-matrix "watermark")

;;
(def watermarked (add fft-matrix-r wm-matrix))

(def new-ifft-matrix (idft-ri-reals watermarked fft-matrix-i))
(display-from-matrix new-ifft-matrix "watermarked")


;; decode here

(def fft-matrix-r-new (dft-reals new-ifft-matrix))

(def decoded-wm (- fft-matrix-r-new fft-matrix-r))
 (img/show (matrix-to-image decoded-wm (fn [x] (* x 1))) :title "decoded")
