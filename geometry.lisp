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
   using formula \rho (x,y) =\sqrt{\Sigma_{i=1}^3 (x_i-y_i)^2}
   Actually, not used in tree processing"
  (declare (optimize (speed 3))
           (type dot dot1 dot2))
  (let ((components (make-array 3
                                :element-type 'single-float
                                :initial-element 0.0)))
    (declare (dynamic-extent components))
    (sqrt
          (reduce #'+ (map-into components
                                #'(lambda (x y)
                                    (declare (type single-float x y))
                                    (expt (- x y) 2))
                                dot1 dot2)))))

(defun hit-box (min max origin dir)
  "Find intersection of a ray and axis-aligned
   cuboid"
  ;; FIXME: Partially copies implementation of dot-betweenp and hit-plane
  (declare (type dot min max origin dir)
           (optimize (speed 3)))
  
  (let ((insidep t)
        (inside (make-array 3))
        (candidate-plane (make-array 3 :element-type 'single-float))
        (tdist (make-array 3 :element-type 'single-float))
        (coord (make-array 3 :element-type 'single-float)))
    (declare (dynamic-extent inside candidate-plane tdist))
    
    (dotimes (i 3)
      (cond
        ((< (aref origin i)
            (aref min i))
         (setf insidep nil
               (aref candidate-plane i) (aref min i)
               (aref inside i) nil))

        ((> (aref origin i)
            (aref max i))
         (setf insidep nil
               (aref candidate-plane i) (aref max i)
               (aref inside i) nil))

        (t
         (setf (aref inside i) t))))

    (if insidep
        (return-from hit-box (values t (copy-seq origin)))) ; COPY-SEQ for safety

    (dotimes (i 3)
      (setf (aref tdist i)
            (if (or (aref inside i)
                    (= (aref dir i) 0.0))
                -1.0
                (/ (- (aref candidate-plane i)
                      (aref origin i))
                   (aref dir i)))))

    (let* ((maxt (reduce #'max tdist))
           (planenum (position maxt tdist)))
      (declare (type single-float maxt))
      (if (< maxt 0.0) (return-from hit-box nil))
      
      (dotimes (i 3)
        (setf (aref coord i)
              (if (= i planenum)
                  (aref candidate-plane i)
                  
                  (let ((coord-i (+ (aref origin i)
                                    (* maxt (aref dir i)))))
                    (if (or (< coord-i (aref min i))
                            (> coord-i (aref max i)))
                        (return-from hit-box nil))
                    coord-i))))
      
      (values t coord))))
        

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
             (dotimes (i 3)
               (setf (aref coord i)
                     (if (= i planenum)
                         (aref planedot planenum)
                         (+ (aref origin i) (* k (aref dir i))))))
             (values t coord)))))))

(defun dot-betweenp (dot min max)
  "T if dot is placed between MIN and MAX"
  (declare (type dot dot min max)
           (optimize (speed 3)))

  (every #'(lambda (x min max)
             (and (>= x min)
                  (<= x max)))
         dot min max))
