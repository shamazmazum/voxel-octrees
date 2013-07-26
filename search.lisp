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

(declaim (type (integer 0 #.most-positive-fixnum) *max-depth*))
(defparameter *max-depth* 2
  "A parameter used by local-ray-tree-intersection.
   The less this parameter is, the less is a penalty
   in the worst case")

;; (Maybe) faster version

(defun construct-subspace-indices (intersects)
  (let ((res (list (first intersects)))
         (cur-idx (caar intersects)))
    (dolist (inter (cdr intersects))
      (setq cur-idx (logxor cur-idx (ash 1 (car inter))))
      (push (cons cur-idx (cdr inter)) res))
    (reverse res)))

(defun ray-tree-intersection (tree origin dir &key save-history history)
  "Search intersection of a ray with closest to the ORIGIN
   cuboid in TREE"
  (declare (optimize (speed 3)))
  (if save-history (setq history (cons tree history)))
  (let* ((bb (node-bounding-box tree))
         (bb-min (car bb))
         (bb-max (cdr bb)))

    (multiple-value-bind (bb-interp bb-inter-coord)
        (and bb (hit-box bb-min bb-max
                         origin dir))
      (if (not bb-interp) (return-from ray-tree-intersection nil))
    
      (cond
        ((leafp tree)
         (let (leaf-intersections
               (dots (node-dots tree)))
           (dotimes (i (length (sb-ext:truly-the simple-vector dots)))
             (multiple-value-bind (interp coord)
                 (let ((dot (aref dots i)))
                   (hit-box dot (sum-vector (copy-seq (sb-ext:truly-the dot dot)) *voxel*)  ;; Extra allocation
                            origin dir))
               (if interp (push coord leaf-intersections))))
           
           (if leaf-intersections
               (flet ((get-closest (dot1 dot2)
                        (if (< (calc-abs-metric dot1 origin)
                               (calc-abs-metric dot2 origin))
                            dot1 dot2)))
                 (list
                  (reduce #'get-closest leaf-intersections)
                  history)))))
        
        (t
         (let ((center (node-dots tree))
               plane-intersections)
           
           (dotimes (i 3)
             (multiple-value-bind (interp coord) (hit-plane origin dir center i)
               (if (and interp (dot-betweenp coord bb-min bb-max))
                   (push (cons i coord) plane-intersections))))

           (setq plane-intersections (construct-subspace-indices
                                      (cons (cons (get-subspace-idx center bb-inter-coord) bb-inter-coord)
                                            (flet ((get-closest-tagged (dot1 dot2)
                                                     (< (calc-abs-metric (cdr dot1) origin)
                                                        (calc-abs-metric (cdr dot2) origin))))
                                              (sort plane-intersections #'get-closest-tagged)))))

           (dolist (subspace plane-intersections)
             (let ((inter (ray-tree-intersection (nth (car subspace)
                                                      (node-children tree))
                                                 (cdr subspace) dir
                                                 :save-history save-history
                                                 :history history)))
               (if inter (return inter))))))))))

;; More precise but slower
;; ?????
#|(defun ray-tree-intersection (tree origin dir)
  "Search intersection of a ray with closest to the ORIGIN
   cuboid in TREE"
  (declare (optimize (speed 3)))

  (cond
    ((leafp tree)
     (let (leaf-intersections
           (dots (node-dots tree)))
       (dotimes (i (length (sb-ext:truly-the simple-vector dots)))
         (multiple-value-bind (interp coord)
             (let ((dot (aref dots i)))
               (hit-box dot (sum-vector (copy-seq (sb-ext:truly-the dot dot)) *voxel*)  ;; Extra allocation
                        origin dir))
           (if interp (push coord leaf-intersections))))

       (if leaf-intersections
           (flet ((get-closest (dot1 dot2)
                    (if (< (calc-abs-metric dot1 origin)
                           (calc-abs-metric dot2 origin))
                        dot1 dot2)))
             (reduce #'get-closest leaf-intersections)))))
        
        (t
         (let ((bb-intersections
                (loop for child in (node-children tree)
                   for bb = (node-bounding-box child)
                   for inter = (if (and bb (hit-box (car bb) (cdr bb) origin dir)) (car bb))
                   when inter collect (cons child inter))))
           
           (setq bb-intersections
                 (flet ((get-closest-tagged (dot1 dot2)
                          (< (calc-abs-metric (cdr dot1) origin)
                             (calc-abs-metric (cdr dot2) origin))))
                   (sort bb-intersections #'get-closest-tagged)))

           (loop for subspace in bb-intersections
              for inter = (ray-tree-intersection (car subspace) origin dir)
              until inter finally (return inter))))))|#



;; The next function is useful for searching intersections
;; with a tree and a group of rays with insignificantly
;; different directions and the same origin.
;; To do this call RAY-TREE-INTERSECTION with the first ray and true :SAVE-HISTORY
;; With next rays call LOCAL-RAY-TREE-INTERSECTION with HISTORY obtained from the call
;; of RAY-TREE-INTERSECTION

;; NB: There is a penalty in speed (with comparsion with RAY-TREE-INTERSECTION),
;; if rays are not local or do not hit the same continuous object

(defun local-ray-tree-intersection (history origin dir &optional (depth *max-depth*))
  "Serves the same purpose as RAY-TREE-INTERSECTION.
   The first argument is a list of nested tree nodes:
   the first node is a child of the second, the second
   is a child of the third and so on. ORIGIN and DIR are
   arrays of type DOT

   Returns intersection and (updated history or NIL if history
   was not updated)"
  (declare (optimize (speed 3))
           (type list history)
           (type (integer 0 #.most-positive-fixnum) depth))
  (if (or (= depth 0)
          (null (cdr history)))
      (ray-tree-intersection (car (last history)) origin dir
                             :save-history t)
      (let ((intersection (ray-tree-intersection (car history) origin dir)))
        (if intersection intersection
            (local-ray-tree-intersection (cdr history) origin dir (sb-ext:truly-the
                                                                   (integer 0 #.most-positive-fixnum)
                                                                   (1- depth)))))))
