module methode_resolution
  implicit none
  !Declaration sur type matrice diagonale
  type matdiag
     real*8,dimension(:,:),allocatable::valeurs
     integer,dimension(:),allocatable::dist_diag
  endtype matdiag

contains
  
  subroutine init(mat) !Routine d'initialisation
    implicit none
    common /noeud/ m !Nombre de noeuds par direction
    common /param/ gb !Alpha et beta
    type(matdiag),intent(inout)::mat
    real*8,dimension(2)::gb
    integer,parameter::nbdiag=5 !La matrice possede uniquement 5 diagonales
    real*8,dimension(nbdiag)::val
    integer::i,j,n,k,m
    real*8::h

    n=m**2 !Nombre de noeuds total
    h=1.D0/(1.D0+m) !Calcul du pas

    allocate(mat%valeurs(n,nbdiag),mat%dist_diag(nbdiag))
    !Affectation des valeurs necessaires pour construire la matrice
    mat%valeurs=0
    mat%dist_diag=(/-m,-1,0,1,m/)
    val=(/-gb(2)*h-2.D0,&
         &-gb(2)*h-2.D0,&
         &8.D0+2*gb(1)*h**2,&
         &gb(2)*h-2.D0,&
         &gb(2)*h-2.D0/)
    !Construction de la matrice
    do i=1,n
       do j=1,nbdiag
          k=i+mat%dist_diag(j)
          if (k.GT.0 .and. k.LE.n) then
             mat%valeurs(i,j)=val(j)/(2*h**2)
          endif
       enddo
    enddo
    !Affectation des zeros qui correspondent au bord
    do i=n-m,m,-m
       mat%valeurs(i,4)=0
    enddo
    do i=m+1,n-m+1,m
       mat%valeurs(i,2)=0
    enddo

  endsubroutine init

  subroutine vectb(b) !Routine pour le second membre
    implicit none
    common /noeud/ m !Nombre de noeuds par direction
    common /param/ gb !Alpha et beta
    real*8,dimension(:),intent(out),allocatable::b
    real*8,dimension(2)::gb
    integer::i,j,m
    real*8::h

    allocate(b(m**2))
    h=1.D0/(1.D0+m) !Calcul du pas
    !b=f(i,j) avec f donnee dans le rapport + reindexation 2D-1D
    do i=1,m
       do j=1,m
          b(i+(j-1)*m)=2*sin(4*acos(0.)*i*h)*exp(-(j*h)**2)*&
               &(1-2*(j*h)**2+8*(acos(0.))**2+&
               &gb(1)/2.D0+&
               &gb(2)*(2*acos(0.)/tan(4*acos(0.)*i*h)-j*h))
       end do
    end do
    !Prise en compte des conditions limites
    do i=1,m
       b(i+(m-1)*m)=b(i+(m-1)*m)+(sin(4*acos(0.)*i*h)*exp(-1.))/h**2
       b(i+(1-1)*m)=b(i+(1-1)*m)+(sin(4*acos(0.)*i*h)*exp(-0.))/h**2
       b(1+(i-1)*m)=b(1+(i-1)*m)+(sin(4*acos(0.)*0)*exp(-(i*h)**2))/h**2
       b(m+(i-1)*m)=b(m+(i-1)*m)+(sin(4*acos(0.)*1)*exp(-(i*h)**2))/h**2
    end do

  endsubroutine vectb

  function norme2(v) !Calcul de la norme 2 d'un vecteur
    implicit none
    real*8,dimension(:)::v
    real*8::norme2
    integer::i

    norme2=0
    do i=1,size(v)
       norme2=norme2+v(i)**2
    enddo

    norme2=sqrt(norme2)

  endfunction norme2

  function prodscal(u,v) !Produit scalaire entre u et v
    implicit none
    real*8,dimension(:)::u,v
    real*8::prodscal
    integer::i

    prodscal=0
    do i=1,size(u)
       prodscal=prodscal+u(i)*v(i)
    enddo

  endfunction prodscal

  subroutine prod(mat,x,b) !Produit de Ax et le stock dans b
    implicit none
    type(matdiag),intent(in)::mat
    real*8,dimension(:),intent(in)::x
    real*8,dimension(:),intent(out),allocatable::b
    real*8,dimension(:),allocatable::x_etendu
    integer::i,j,n,m

    n=size(mat%valeurs(:,1))
    m=size(mat%valeurs(1,:))
    !On a besoin d'un vecteur x_etendu pour eviter de multiplier par une
    !valeur hors tableau (a cause du stockage diagonal)
    allocate(x_etendu(mat%dist_diag(1)+1:n+mat%dist_diag(m)),b(n))
    x_etendu=0
    x_etendu(1:n)=x
    b=0
    !Produit matrice-vecteur dans le cas d'un stockage diagonal
    do i=1,n
       do j=1,m     
          b(i)=b(i)+mat%valeurs(i,j)*x_etendu(i+mat%dist_diag(j))
       enddo
    enddo

  endsubroutine prod

  subroutine prodtrans(mat,x,b) !Produit de Atx et le stock dans b
    implicit none
    type(matdiag),intent(in)::mat
    real*8,dimension(:),intent(in)::x
    real*8,dimension(:),intent(out),allocatable::b
    real*8,dimension(:),allocatable::x_etendu
    integer::i,j,n,m

    n=size(mat%valeurs(:,1))
    m=size(mat%valeurs(1,:))
    !Meme raison pour le vecteur etendu
    allocate(x_etendu(mat%dist_diag(1)+1:n+mat%dist_diag(m)),b(n))
    x_etendu=0
    x_etendu(1:n)=x
    b=0
    !Produit matrice transposee-vecteur dans le cas d'un stockage diagonal,
    !sans avoir la matrice transposee de stockee
    do i=1,n
       b(i)=b(i)+mat%valeurs(i,3)*x_etendu(i)
       do j=1,(m-1)/2
          b(i)=b(i)+mat%valeurs(i+mat%dist_diag(j),m+1-j)*x_etendu(i+mat%dist_diag(j))
          b(i)=b(i)+mat%valeurs(i+mat%dist_diag(m+1-j),j)*x_etendu(i+mat%dist_diag(m+1-j))
       enddo
    enddo

  endsubroutine prodtrans

  subroutine grad(mat,x,b) !Gradient simple
    implicit none
    common /export/ out
    common /file/ nom
    type(matdiag),intent(in)::mat
    real*8,dimension(:),intent(in)::b
    real*8,dimension(:),intent(out),allocatable::x
    real*8,dimension(:),allocatable::r,Ar
    real*8::alp,epsi,normerb,normeb,normer,t1,t2
    integer::n,cpt,i
    character(len=30)::nom
    logical::out

    if (out) then !Si le parametre export est vrai alors on cree le fichier
       open(11,file=nom)
    endif

    print*,"Solution approchee par methode du gradient simple : "
    n=size(b)
    allocate(r(n),x(n),Ar(n))
    epsi=epsilon(epsi) !Fonction intrinseque de Fortran pour le calcul d'epsi
    print*,"L'epsilon machine vaut : ",epsi
    cpt=0
    x=0
    normerb=1 !Initialisation pour rentrer au moins une fois dans la boucle
    normeb=norme2(b) !On stock la norme de b car ne change pas au cours du temps

    call cpu_time(t1)
    !voir algo du rapport
    call prod(mat,x,r)
    do i=1,n
       r(i)=b(i)-r(i)
    enddo
    do while (normerb .GE. epsi)
       call prod(mat,r,Ar)
       alp=prodscal(r,r)/prodscal(Ar,r)
       do i=1,n
          x(i)=x(i)+alp*r(i)
          r(i)=r(i)-alp*Ar(i)
       enddo
       normer=norme2(r)
       normerb=normer/normeb
       cpt=cpt+1
       if (out) then !Si le parametre export est vrai alors on ecrit
          write(11,'(I15,3D30.16)'),cpt,normerb,erreur_absolue(x)
       endif
    enddo
    call cpu_time(t2)
    !Information principale liee aux resultats
    print*,"Nombre d'iterations : ",cpt
    print*,"L'erreur relative est de : ",erreur_relative(x)
    print*,"L'erreur absolue est de : ",erreur_absolue(x)
    print*,"Solution calculee en : ",t2-t1,"secondes."

    if (out) then !Si le parametre export est vrai, on ferme le fichier
       close(11)
    endif

  endsubroutine grad

  subroutine gradconju(mat,x,b) !Gradient conjugue, meme commentaire que le gradient simple + voir algo
    implicit none
    common /export/ out
    common /file/ nom
    type(matdiag),intent(in)::mat
    real*8,dimension(:),intent(in)::b
    real*8,dimension(:),intent(out),allocatable::x
    real*8,dimension(:),allocatable::r,Ap,p,rk
    real*8::alp,epsi,normerb,normeb,normer,beta,t1,t2
    integer::n,cpt,i
    character(len=30)::nom
    logical::out

    if (out) then
       open(11,file=nom)
    endif
    
    print*,"Solution approchee par methode du gradient conjugue : "
    n=size(b)
    allocate(r(n),rk(n),x(n),p(n),Ap(n))
    epsi=epsilon(epsi)
    print*,"L'epsilon machine vaut : ",epsi
    cpt=0
    x=0
    normerb=1
    normeb=norme2(b)

    call cpu_time(t1)
    call prod(mat,x,r)
    do i=1,n
       r(i)=b(i)-r(i)
       p(i)=r(i)
    enddo
    do while (normerb .GE. epsi)
       call prod(mat,p,Ap)
       alp=prodscal(r,r)/prodscal(Ap,p)
       do i=1,n
          x(i)=x(i)+alp*p(i)
          rk(i)=r(i)-alp*Ap(i)
       enddo
       beta=prodscal(rk,rk)/prodscal(r,r)
       do i=1,n
          p(i)=rk(i)+beta*p(i)
       enddo
       r=rk
       normer=norme2(r)
       normerb=normer/normeb
       cpt=cpt+1
       if (out) then
          write(11,'(I15,3D30.16)'),cpt,normerb,erreur_absolue(x)
       endif
    enddo
    call cpu_time(t2)
    print*,"Nombre d'iterations : ",cpt
    print*,"L'erreur relative est de : ",erreur_relative(x)
    print*,"L'erreur absolue est de : ",erreur_absolue(x)
    print*,"Solution calculee en : ",t2-t1,"secondes."

    if (out) then
       close(11)
    endif

  endsubroutine gradconju

  subroutine gmres(mat,x,b) !GmRes, meme commentaire que le gradient simple + voir algo
    implicit none
    common /export/ out
    common /file/ nom
    type(matdiag),intent(in)::mat
    real*8,dimension(:),intent(in)::b
    real*8,dimension(:),intent(out),allocatable::x
    real*8,dimension(:),allocatable::r,p,rk,q,Ar
    real*8::alp,epsi,normerb,normeb,normer,beta,t1,t2
    integer::n,cpt,i
    character(len=30)::nom
    logical::out

    if (out) then
       open(11,file=nom)
    endif

    print*,"Solution approchee par methode du residu minimal : "
    n=size(b)
    allocate(r(n),rk(n),x(n),p(n),q(n),Ar(n))
    epsi=epsilon(epsi)
    print*,"L'epsilon machine vaut : ",epsi
    cpt=0
    x=0
    normerb=1
    normeb=norme2(b)

    call cpu_time(t1)
    call prod(mat,x,r)
    do i=1,n
       r(i)=b(i)-r(i)
       p(i)=r(i)
    enddo
    call prod(mat,p,q)
    do while (normerb .GE. epsi)
       alp=prodscal(r,q)/prodscal(q,q)
       do i=1,n
          x(i)=x(i)+alp*p(i)
          rk(i)=r(i)-alp*q(i)
       enddo
       call prod(mat,rk,Ar)
       beta=-prodscal(Ar,q)/prodscal(q,q)
       do i=1,n
          p(i)=rk(i)+beta*p(i)
          q(i)=Ar(i)+beta*q(i)
       enddo
       r=rk
       normer=norme2(r)
       normerb=normer/normeb
       cpt=cpt+1
       if (out) then
          write(11,'(I15,3D30.16)'),cpt,normerb,erreur_absolue(x)
       endif
    enddo
    call cpu_time(t2)
    print*,"Nombre d'iterations : ",cpt
    print*,"L'erreur relative est de : ",erreur_relative(x)
    print*,"L'erreur absolue est de : ",erreur_absolue(x)
    print*,"Solution calculee en : ",t2-t1,"secondes."

    if (out) then
       close(11)
    endif

  endsubroutine gmres

  subroutine bicg(mat,x,b) !BiCG, meme commentaire que le gradient simple + voir algo
    implicit none
    common /export/ out
    common /file/ nom
    type(matdiag),intent(in)::mat
    real*8,dimension(:),intent(in)::b
    real*8,dimension(:),intent(out),allocatable::x
    real*8,dimension(:),allocatable::r,rtilde,p,ptilde,rpre,rpretilde,xtilde,btilde,Ap,Atptilde
    real*8::alp,epsi,normerb,normeb,beta,t1,t2
    integer::n,cpt,i
    character(len=30)::nom
    logical:: out

    if (out) then
       open(11,file=nom)
    endif

    print*,"Solution approchee par methode du gradient biconjugue : "
    n=size(b)
    allocate(r(n),rtilde(n),x(n),p(n),ptilde(n),rpre(n),rpretilde(n),xtilde(n),btilde(n),Ap(n),Atptilde(n))
    epsi=epsilon(epsi)
    print*,"L'epsilon machine vaut : ",epsi
    cpt=0
    x=0
    xtilde=0
    btilde=b
    normerb=1
    normeb=norme2(b)

    call cpu_time(t1)
    call prod(mat,x,rpre)
    call prod(mat,xtilde,rpretilde)
    do i=1,n
       rpre(i)=b(i)-rpre(i)
       p(i)=rpre(i)
       rpretilde(i)=btilde(i)-rpretilde(i)
       ptilde(i)=rpretilde(i)
    enddo
    do while (normerb .GE. epsi)
       call prod(mat,p,Ap)
       call prodtrans(mat,ptilde,Atptilde)
       alp=prodscal(rpretilde,rpre)/prodscal(ptilde,Ap)
       do i=1,n
          x(i)=x(i)+alp*p(i)
          r(i)=rpre(i)-alp*Ap(i)
          rtilde(i)=rpretilde(i)-alp*Atptilde(i)
       enddo
       beta=prodscal(rtilde,r)/prodscal(rpretilde,rpre)
       do i=1,n
          p(i)=r(i)+beta*p(i)
          ptilde(i)=rtilde(i)+beta*ptilde(i)
          rpre(i)=r(i)
          rpretilde(i)=rtilde(i)
       enddo
       normerb=norme2(r)/normeb
       cpt=cpt+1
       if (out) then
          write(11,'(I15,3D30.16)'),cpt,normerb,erreur_absolue(x)
       endif
    enddo
    call cpu_time(t2)
    print*,"Nombre d'iterations : ",cpt
    print*,"L'erreur relative est de : ",erreur_relative(x)
    print*,"L'erreur absolue est de : ",erreur_absolue(x)
    print*,"Solution calculee en : ",t2-t1,"secondes."

    if (out) then
       close(11)
    endif

  endsubroutine bicg

  function erreur_relative(x) !Calcul de l'erreur relative d'un vecteur x en norme L1
    implicit none
    common /noeud/ m
    real*8,dimension(:)::x
    real*8::erreur_relative,haut,bas,h,u
    integer::i,j,m

    h=1./(1.+m)
    haut=0
    bas=0
    do i=1,m
       do j=1,m
          u=sin(4*acos(0.)*i*h)*exp(-(j*h)**2)
          haut=haut+abs(x(i+(j-1)*m)-u)
          bas=bas+abs(u)
       enddo
    enddo
    
    erreur_relative=(haut/bas)*h**2 !On multiplie par h**2 car nous sommes en 2D

  endfunction erreur_relative

  function erreur_absolue(x) !Calcul de l'erreur absolue d'un vecteur x en norme L1
    implicit none
    common /noeud/ m
    real*8,dimension(:)::x
    real*8::erreur_absolue,h,u
    integer::i,j,m

    h=1./(1.+m)
    erreur_absolue=0
    do i=1,m
       do j=1,m
          u=sin(4*acos(0.)*i*h)*exp(-(j*h)**2)
          erreur_absolue=erreur_absolue+abs(x(i+(j-1)*m)-u)
       enddo
    enddo

    erreur_absolue=erreur_absolue*h**2 !On multiplie par h**2 car nous sommes en 2D

  endfunction erreur_absolue

endmodule methode_resolution
