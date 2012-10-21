(require :gsll)
(require :antik)
(in-package :antik-user)

;; (defun vequalp (v1 v2)
;;   (iter (for e1 :vector-element v1)
;;         (for e2 :vector-element v2)
;;         (when (/= e1 e2) (leave nil))
;;         (finally (return t))))

(deftype vec-double () '(simple-array double-float (3)))

(defun load-data (&optional (fname "/home/mcstar/pftpfc/gamma.original"))
  (let ((data (with-open-file (s fname)
		(loop for line = (read-line s nil nil)
		   while line collect (read-from-string (format nil "(~a)" line))))))
    (loop for row in data
          collect (make-array 3 :element-type '(signed-byte 8) :initial-contents (sort (butlast row) #'<))
            into directions
          collect (first (last row))
            into gammas
          finally (return (values directions gammas)))))

(defun gen-all-directions (dirs)
  (let ((all-dirs (make-hash-table :test #'equalp))
        (indices (make-simple-grid :grid-type 'foreign-array
                                   :element-type '(signed-byte 8)
                                   :initial-contents '(0 1 2)))
        (p (make-permutation 3))
        (signs (list #(1 1 1) #(-1 1 1) #(1 -1 1) #(1 1 -1) #(-1 -1 1) #(-1 1 -1) #(1 -1 -1) #(-1 -1 -1))))
    (let ((pdirs (loop 
                   for i0 = (aref indices 0)
                   for i1 = (aref indices 1)
                   for i2 = (aref indices 2)
                   append (loop for d in dirs collect
                     (make-array 3
                          :element-type '(signed-byte 8)
                          :initial-contents (list (aref d i0)
                                                  (aref d i1)
                                                  (aref d i2))))
                   while (multiple-value-bind (_ end)
                             (permutation-next p)
                           (declare (ignore _)) end)
                   do (permute p indices))))
      (loop for sign in signs
            do (loop for d in pdirs
                     do (setf (gethash (* sign d) all-dirs) t))))
    all-dirs))
    
(defun getGamma (vec)
  (let ((dir (sort (abs vec) #'<)))
    (gethash dir *ght*)))

(defun less-than-piover2? (v1 v2)
  (< (realpart
      (acos
       (/ (inner v1 v2)
          (sqrt( * (inner v1 v1)
                   (inner v2 v2))))))
     (/ pi 2.0)))

(defun is-outside? (tv rv)
  (let ((diff (- tv rv)))
    (let ((a2 (inner tv tv))
	  (b2 (inner diff diff)))
      (> (realpart
          (acos
           (/ (+ a2 b2 (- (inner rv rv)))
              (* 2.0 (sqrt (* a2 b2))))))
         (+ (/ pi 2.0) 0.001)))))

(defun gentriples (dirs)
  (let ((boolmat (make-array (list (length dirs) (length dirs))
                             :element-type boole
                             :initial-element nil)))
    
  (let ((len (length vec-list))
	(vecs (apply 'vector vec-list)))
    (loop for i from 0 below len
       for i-el = (aref vecs i)
       append (loop for j from (1+ i) below len
		 for j-el = (aref vecs j)
		 append (loop for k from (1+ j) below len
			   for k-el = (aref vecs k)
			   when (and (less-than-piover2? i-el j-el)
				     (less-than-piover2? i-el k-el)
				     (less-than-piover2? j-el k-el))
			   collect (vector i-el j-el k-el))))))

(defun solve (mat)
  (let ((row0 (aref mat 0))
	(row1 (aref mat 1))
	(row2 (aref mat 2)))
    (let ((gvec (grid:grid (accessgamma row0)
			   (accessgamma row1)
			   (accessgamma row2))))
      (let ((invnorm0 (/ 1.0 (sqrt (grid:inner row0 row0))))
	    (invnorm1 (/ 1.0 (sqrt (grid:inner row1 row1))))
	    (invnorm2 (/ 1.0 (sqrt (grid:inner row2 row2)))))
	(let ((row0 (times-s-v invnorm0 row0))
	      (row1 (times-s-v invnorm1 row1))
	      (row2 (times-s-v invnorm2 row2)))
	  (let ((mat (grid:transpose (grid:matrix-from-columns row0 row1 row2))))
	    (multiple-value-bind (invmat det)
		(handler-case (antik:invert-matrix mat)
		  (error () (values nil nil)))
	      (when det
		(if (< (abs det) 1e-8)
		  nil
		  (let ((tmp (grid:grid 0.0 0.0 0.0)))
		    (vector (inner-v-v (grid:row invmat 0 tmp) gvec)
			    (inner-v-v (grid:row invmat 1 tmp) gvec)
			    (inner-v-v (grid:row invmat 2 tmp) gvec))))))))))))
      

(defun solve-intersection (triple-vecs)
  (loop for one-triple in triple-vecs
       for solution = (solve one-triple)
       when solution
       collect solution))


(defun erase-points (gammas points)
  (loop for point in points
     do
       (setf point (coerce point 'vec-double))
     when (loop for gamma in gammas
	     do
	       (setf gamma (coerce gamma 'vec-double))
	       (when (less-than-piover2? gamma point)
		 (when (or (> (sqrt (inner-v-v point point)) 5.0) (is-outside? gamma point))
		   (return nil)))
	     finally (return t))
     collect point))


(multiple-value-bind (data gammas) (load-data)
  (defparameter *data* data)
  (defparameter *gammas* gammas))

(defparameter *ght* 
  (let ((ht (make-hash-table :test 'equalp)))
    (loop for dir in *data*
          for gamma in *gammas*
          do (setf (gethash dir ht) gamma))
    ht))

(defparameter *all-dirs* (gen-all-directions *data*))

(defparameter *triplets* (gentriples *alldirs*))

(defparameter *points* (solve-intersection *triplets*))

(defparameter *gammas* (loop for dir in *alldirs*
			  for norm = (sqrt (inner-v-v dir dir))
			  collect (times-s-v (* (accessgamma dir) (/ 1.0 norm)) dir)))
     
(defparameter *result* (erase-points *gammas* *points*))


(defun save-for-qhull (points &optional (fname "/home/mcstar/pftpfc/lisp-result"))
  (with-open-file (s fname :direction :output :if-exists :supersede :if-does-not-exist :create)
    (format s "3~%~d~%" (length points))
    (loop for point in points
	 do (format s "~{~a~t~}~%" (coerce point 'list)))))

