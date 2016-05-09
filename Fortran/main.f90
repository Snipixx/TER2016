program main
  use methode_resolution
  implicit none
  character(len=30)::arg1,arg2,arg3,nom
  type(matdiag)::A
  real*8,dimension(:),allocatable::b,x
  real*8,dimension(2)::gb
  integer::m,e
  logical::isdigit,out
  !Variable globale utilisee dans le module
  common /noeud/ m
  common /export/ out
  common /file/ nom
  common /param/ gb
  !Nombre d'argument entre correct
  if ((iargc().GE.1) .AND. (iargc().LE.3)) then
     call getarg(1,arg1) !On recupere les arguments entres
     call getarg(2,arg2)
     call getarg(3,arg3)
     out=(arg3=="-export").OR.(arg3=="-o") !Exportation ou non
     read(arg2,*,IOSTAT=e)m !Lecture du nombre de noeuds par direction
     isdigit=e==0 !Verification si c'est un entier
     if (arg1.NE."-help") then !On demande alpha et beta
        print*,"Entrer la valeur du parametre Beta :"
        read*,gb(2)
        print*,"Entrer la valeur du parametre Gamma :"
        read*,gb(1)
     endif
     nom=trim(arg1(2:30))//"_"//trim(arg2)//".dat" !Nom du ficher d'export
     !Si tout est OK, on applique la methode demandee
     if ((arg1 == "-gs") .AND. isdigit) then
        call init(A)
        call vectb(b)
        call grad(A,x,b)
        deallocate(A%valeurs,A%dist_diag,x,b)
     elseif ((arg1 == "-gc") .AND. isdigit) then
        call init(A)
        call vectb(b)
        call gradconju(A,x,b)
        deallocate(A%valeurs,A%dist_diag,x,b)
     elseif ((arg1 == "-gmres") .AND. isdigit) then
        call init(A)
        call vectb(b)
        call gmres(A,x,b)
        deallocate(A%valeurs,A%dist_diag,x,b)
     elseif ((arg1 == "-bicg") .AND. isdigit) then
        call init(A)
        call vectb(b)
        call bicg(A,x,b)
        deallocate(A%valeurs,A%dist_diag,x,b)
     elseif (arg1 == "-help") then !Affichage de l'help si demandee
        print*,"Liste des commandes disponible :"
        print*,"-gs X // Methode du Gradient simple avec X noeuds par direction"
        print*,"-gc X // Methode du Gradient conjugue avec X noeuds par direction"
        print*,"-gmres X // Methode du Residu minimal avec X noeuds par direction"
        print*,"-bicg X // Methode du Gradient biconjugue avec X noeuds par direction"
        print*,"-export ou -o // Exporter les donnees dans un fichier (a mettre en dernier parametre)"
        print*,"Exemple d'input correct : ./solveur_gfortran -bicg 300"
     else
        print*,"Mauvais argument entre, -help pour la liste des commandes."
     endif
  else !Si trop d'arguments entres
     print*,"Mauvais nombre d'arguments entres, -help pour la liste des commandes."
  endif

endprogram main
