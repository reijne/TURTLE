      Subroutine upCase(line)
      Implicit Integer (a-z)
      Character*(*) line
      Do 100 k=1,len(line)
         Index=iChar(line(k:k))
         If(97.le.Index .and. Index.le.122) line(k:k)=Char(Index-32)
100   Continue
      Return
      End
