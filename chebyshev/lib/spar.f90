module sparsekit
implicit none

contains
!------------------------------------------------------------------------
subroutine coicsr (n,nnz,job,a,ja,ia,iwk)
integer n,nnz,job,i,j,k,inext,jnext,init,ipos
integer ia(nnz),ja(nnz),iwk(n+1)
complex(8) a(*)
!------------------------------------------------------------------------
!IN-PLACE coo-csr conversion routine.
!------------------------------------------------------------------------
!this subroutine converts a matrix stored in coordinate format into
!the csr format. The conversion is done in place in that the arrays
!a,ja,ia of the result are overwritten onto the original arrays.
!------------------------------------------------------------------------
!on entry:
!---------
!n	= integer. row dimension of A.
!nnz	= integer. number of nonzero elements in A.
!job   = integer. Job indicator. when job=1, the real values in a are
!        filled. Otherwise a is not touched and the structure of the
!        array only (i.e. ja, ia)  is obtained.
!a	= real array of size nnz (number of nonzero elements in A)
!        containing the nonzero elements
!ja	= integer array of length nnz containing the column positions
!	  of the corresponding elements in a.
!ia	= integer array of length nnz containing the row positions
!	  of the corresponding elements in a.
!iwk	= integer work array of length n+1
!on return:
!----------
!a
!ja
!ia	= contains the compressed sparse row data structure for the
!        resulting matrix.
!Note:
!-------
!        the entries of the output matrix are not sorted (the column
!        indices in each are not in increasing order) use coocsr
!        if you want them sorted.
!----------------------------------------------------------------------c
! Coded by Y. Saad, Sep. 26 1989                                      c
!----------------------------------------------------------------------c
complex(8) t,tnext
logical values
!-----------------------------------------------------------------------
values = (job .eq. 1)
!find pointer array for resulting matrix.
do 35 i=1,n+1
iwk(i) = 0
35   continue
do 4 k=1,nnz
i = ia(k)
iwk(i+1) = iwk(i+1)+1
4    continue
!------------------------------------------------------------------------
iwk(1) = 1
do 44 i=2,n
iwk(i) = iwk(i-1) + iwk(i)
44   continue
!
!    loop for a cycle in chasing process.
!
init = 1
k = 0
5    if (values) t = a(init)
i = ia(init)
j = ja(init)
ia(init) = -1
!------------------------------------------------------------------------
6    k = k+1
!    current row number is i.  determine  where to go.
ipos = iwk(i)
!    save the chased element.
if (values) tnext = a(ipos)
inext = ia(ipos)
jnext = ja(ipos)
!    then occupy its location.
if (values) a(ipos)  = t
ja(ipos) = j
!    update pointer information for next element to come in row i.
iwk(i) = ipos+1
!    determine  next element to be chased,
if (ia(ipos) .lt. 0) goto 65
t = tnext
i = inext
j = jnext
ia(ipos) = -1
if (k .lt. nnz) goto 6
goto 70
65   init = init+1
if (init .gt. nnz) goto 70
if (ia(init) .lt. 0) goto 65
!    restart chasing --
goto 5
70   do 80 i=1,n
ia(i+1) = iwk(i)
80   continue
ia(1) = 1
return
!----------------- end of coicsr ----------------------------------------
!------------------------------------------------------------------------
end

      subroutine amux (n, x, y, a,ja,ia) 
      complex(8)  x(*), y(*), a(*)
      integer n, ja(*), ia(*)
!-----------------------------------------------------------------------
!        A times a vector
!----------------------------------------------------------------------- 
!multiplies a matrix by a vector using the dot product form
!Matrix A is stored in compressed sparse row storage.
!
!on entry:
!----------
!n     = row dimension of A
!x     = real array of length equal to the column dimension of
!        the A matrix.
!a, ja,
!   ia = input matrix in compressed sparse row format.
!
!on return:
!-----------
!y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
!local variables
!
      complex(8) t
      integer i, k
!-----------------------------------------------------------------------
      do 100 i = 1,n
!
!    compute the inner product of row i with vector x
!
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
 99      continue
!
!    store result in y(i) 
!
         y(i) = t
 100  continue
!
      return
!---------end-of-amux---------------------------------------------------
!-----------------------------------------------------------------------
      end

      subroutine atmux (n, x, y, a, ja, ia)
      complex(8) x(*), y(*), a(*)
      integer n, ia(*), ja(*)
!-----------------------------------------------------------------------
!        transp( A ) times a vector
!----------------------------------------------------------------------- 
!multiplies the transpose of a matrix by a vector when the original
!matrix is stored in compressed sparse row storage. Can also be
!viewed as the product of a matrix by a vector when the original
!matrix is stored in the compressed sparse column format.
!-----------------------------------------------------------------------
!
!on entry:
!----------
!n     = row dimension of A
!x     = real array of length equal to the column dimension of
!        the A matrix.
!a, ja,
!   ia = input matrix in compressed sparse row format.
!
!on return:
!-----------
!y     = real array of length n, containing the product y=transp(A)*x
!
!-----------------------------------------------------------------------
!    local variables 
!
      integer i, k 
!-----------------------------------------------------------------------
!
!    zero out output vector
!
      do 1 i=1,n
         y(i) = 0.0
 1    continue
!
!loop over the rows
!
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue
!
      return
!-------------end-of-atmux---------------------------------------------- 
!-----------------------------------------------------------------------
      end
end module sparsekit

