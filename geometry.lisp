;; Copyright (c) 2013, Vasily Postnicov
;; All rights reserved.
;;
;; Redistribution and use in source and binary forms, with or without modification,
;; are permitted provided that the following conditions are met:

;;  Redistributions of source code must retain the above copyright notice, this
;;  list of conditions and the following disclaimer.
;;
;;  Redistributions in binary form must reproduce the above copyright notice, this
;;  list of conditions and the following disclaimer in the documentation and/or
;;  other materials provided with the distribution.
;;
;;THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
;;ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
;;WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
;;DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
;;ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
;;(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
;;LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
;;ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

(in-package :voxel-octrees)

#|(defun calc-dist (dot1 dot2)
  (reduce #'+ (map 'vector #'(lambda (x y)
                               (abs (- x y)))
                   dot1 dot2)))|#

(declaim (ftype (function (dot dot) single-float) calc-abs-metric))
(defun calc-abs-metric (dot1 dot2)
  "Calculate distance between dot1 and dot2
   using formula \rho(x,y)= \Sigma_{i=1}^3 \lvert x_i - y_i \rvert"
  (declare (optimize (speed 3))
           (type dot dot1 dot2))
  (let ((components (make-array 3
                                :element-type 'single-float
                                :initial-element 0.0)))
    (declare (dynamic-extent components))
    (reduce #'+ (map-into components
                          #'(lambda (x y)
                              (declare (type single-float x y))
                              (abs (- x y)))
                          dot1 dot2))))

(declaim (ftype (function (dot dot) single-float) calc-sqr-metric))
(defun calc-sqr-metric (dot1 dot2)
  "Calculate distance between dot1 and dot2
   using formula \rho (x,y) = \Sigma_{i=1}^3 (x_i-y_i)^2

   A square of usual euclid metric."
  (declare (optimize (speed 3))
           (type dot dot1 dot2))
  (let ((components (make-array 3
                                :element-type 'single-float
                                :initial-element 0.0)))
    (declare (dynamic-extent components))
    (reduce #'+ (map-into components
                          #'(lambda (x y)
                              (declare (type single-float x y))
                              (expt (- x y) 2))
                          dot1 dot2))))

(declaim (ftype (function (dot dot dot &optional dot) (values dot boolean))
                fit-into-box))
(defun fit-into-box (min max dot &optional (res (copy-seq dot)))
  "Returns a dot laying in (or on the surface of) the
   axis-aligned box defined by the coordinates MIN
   and MAX and which is the closest dot to DOT.
   Also returns the second value which is T if the
   first returning value and DOT is the same"
  (declare (type dot min max dot res)
           (optimize (speed 3)))
  
  (flet ((choose-closest (min max coord)
           (declare (type single-float min max coord))
           (cond
             ((< coord min) min)
             ((> coord max) max)
             (t coord))))

    (values
     (map-into res #'choose-closest min max dot)
     (equalp res dot))))

(defun dot-betweenp (dot min max)
  "T if dot is placed between MIN and MAX"
  (declare (type dot dot min max)
           (optimize (speed 3)))

  (every #'(lambda (x min max)
             (and (>= x min)
                  (<= x max)))
         dot min max))

;; The following function was taken from C Graphics Gems
;; Of course, it would be better, if it (was) implemented in C

(defun hit-box (min max origin dir)
  "Find intersection of a ray and axis-aligned
   cuboid"
  (declare (optimize (speed 3)))
  
  (let ((tdist (make-array 3 :element-type 'single-float))
        (candidate-plane (make-array 3 :element-type 'single-float))
        (coord (make-array 3 :element-type 'single-float)))
    (declare (dynamic-extent tdist candidate-plane)
             (type dot origin dir tdist coord))

    (multiple-value-bind (candidate-plane insidep)
        (fit-into-box min max origin candidate-plane)

      (cond
        (insidep (values t (copy-seq origin)))
        (t
         (flet ((calc-tdist (origin candidate dir)
                  (if (or (= candidate origin)
                          (= dir 0.0))
                      -1.0
                      (/ (- candidate origin) dir))))
           
           (map-into tdist #'calc-tdist origin candidate-plane dir))
         
         (let* ((maxt (reduce #'max tdist))
                (planenum (position maxt tdist)))
           (declare (type single-float maxt))
           
           (if (>= maxt 0.0)
               (progn
                 (flet ((calc-hit-coords (candidate origin dir i)
                          (declare (type fixnum i))
                          (if (= i planenum) candidate
                              (+ origin (* maxt dir)))))
                   
                   (map-into coord #'calc-hit-coords candidate-plane origin dir '(0 1 2)))

                 (if (dot-betweenp coord min max) (values t coord))))))))))
      
(defun hit-plane (origin dir planedot planenum)
  "Find intersection of a ray and a plane
   with constant coordinate PLANENUM
   and the dot PLANEDOT laying on that plane"
  (declare (type dot origin dir planedot)
           (type fixnum planenum)
           (optimize (speed 3)))
  (cond
    ((= (aref dir planenum) 0.0)
     ;; XXX: Ignore the case of equality of origin[planenum]
     ;; and planedot[planenum] because this does not flip space
     ;; index bit
     nil)
    
    (t
     (let ((k (/ (- (aref planedot planenum)
                    (aref origin planenum))
                 (aref dir planenum))))
       (if (< k 0.0) nil ; Case k = 0 is above
           (let ((coord (make-array 3 :element-type 'single-float)))
             (flet ((calc-hit-coords (planecoord origin dir i)
                      (declare (type fixnum i))
                      (if (= i planenum) planecoord
                          (+ origin (* k dir)))))
               (values t (map-into coord #'calc-hit-coords planedot origin dir '(0 1 2))))))))))

(declaim (ftype (function (dot dot dot single-float) boolean) box-ball-interp))
(defun box-ball-interp (min max center radius)
  "Checks if a axis-aligned box hits a ball."
  (declare (optimize (speed 3))
           (type single-float radius))

  (let ((fit-center (make-array 3 :element-type 'single-float)))
    (declare (dynamic-extent fit-center))
    (< (calc-sqr-metric
        (fit-into-box min max center fit-center)
        center)
       (expt radius 2))))

(defun closest-in-set (dot set metric)
  (declare (type dot dot)
           (type (function (dot dot) single-float) metric)
           (optimize (speed 3)))
  
  (flet ((closest-between-two (d1 d2)
           (if (< (funcall metric dot d1)
                  (funcall metric dot d2))
               d1 d2)))
    (reduce #'closest-between-two set)))
