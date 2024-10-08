MODULE save_file
!Insieme di funzioni e subroutine necessarie per il salvataggio dei dati e la scrittura su file
IMPLICIT NONE

    CONTAINS

    INTEGER FUNCTION rows_count(filename) !Si occupa del conteggio delle righe di un file
        CHARACTER(LEN=*), INTENT(IN) :: filename

        rows_count=0
        OPEN(12,file = filename)
        DO
            READ(12,*,end=100)
            rows_count=rows_count+1
        END DO
        100 CONTINUE
        CLOSE(12)

        !Sottrae le righe dove non sono presenti dati
        IF(filename == 'newtable2.txt') THEN
            rows_count = rows_count - 31
            PRINT *, "Sono state lette ", rows_count, " righe dal file '"//filename//"'."
            PRINT *, ""
        ELSE
            rows_count = rows_count - 5        
        END IF
    END FUNCTION rows_count

    LOGICAL FUNCTION monotone_check(vect,dim) !Dato un vettore di magnitudini restituisce se la curva di luce è monotona
        INTEGER, INTENT(IN) :: dim
        REAL*8, DIMENSION(dim), INTENT(IN) :: vect

        INTEGER :: i

        monotone_check = .false. !La supernova viene assunta monotona
        DO i = 1,SIZE(vect)
            IF (vect(i).le.99.0d0) THEN !Cerca il primo valore effettivamente misurato
                IF (vect(i) /= MINVAL(vect)) THEN
                    monotone_check = .true. !Se la condizione viene soddisfatta la supernova non è monotona
                END IF
                EXIT
            END IF
        END DO
    END FUNCTION monotone_check

    SUBROUTINE newtable_buildmat(newtable_name, n_sn, newtable, names) !Si occupa di salvare i dati richiesti del file 'newtable2.txt' nelle matrici 'new_table' e 'names'
        CHARACTER(LEN=*), INTENT(IN) :: newtable_name
        INTEGER, INTENT(OUT) :: n_sn
        REAL*8, ALLOCATABLE, INTENT(OUT) :: newtable(:,:)
        CHARACTER(LEN=6), ALLOCATABLE, INTENT(OUT) :: names(:)

        INTEGER :: i, j
        REAL*8 :: temp

        n_sn = rows_count(newtable_name)

        ALLOCATE(names(n_sn),newtable(n_sn,3))

        OPEN(10,file=newtable_name)
        DO i = 1,31 !Salta le prime 31 righe
            READ(10,*)
        END DO

        DO i=1,n_sn !Inizia la lettura dei dati nel file
            READ(10,*) names(i), temp, newtable(i,1), (temp, j = 1,6), newtable(i,2), temp, newtable(i,3)
        END DO
        CLOSE(10)
    END SUBROUTINE newtable_buildmat

    SUBROUTINE sn_buildmat(n_sn, names, sn_mat) !Si occupa di salvare i dati richiesti delle singole supernovae nella matrice 'sn_mat'
        INTEGER, INTENT(IN) :: n_sn
        CHARACTER(LEN=6), INTENT(IN) :: names(n_sn)
        REAL*8, ALLOCATABLE, INTENT(OUT) :: sn_mat(:,:)
        
        INTEGER, DIMENSION(n_sn) :: len_arr
        INTEGER :: max_len, posB
        INTEGER :: i, j, k, n, t
        REAL*8 :: tempMJD, tempB, tempBerr
        CHARACTER(LEN=1) :: str

        max_len = 0
        DO i = 1,n_sn !Sceglie di quante righe allocare la matrice 'sn_mat'
            len_arr(i) = rows_count("SN"//names(i)//".dat")
            IF(len_arr(i).gt.max_len) max_len = len_arr(i)
        END DO

        ALLOCATE(sn_mat(max_len,n_sn*3))
        sn_mat = 99.9d0
        k = 1
        PRINT *, "Salvataggio dei dati in corso..."
        DO n = 1,n_sn
            j = 1
            posB = 3
            OPEN(11,file = "SN"//names(n)//".dat")
            DO i = 1, len_arr(n)+5
                IF(i.lt.5) THEN !Salta le prime 4 righe
                    READ(11,*)
                ELSE IF(i.eq.5) THEN !Controlla qual'è la colonna della magnitudine 'B'
                    READ(11,*) (str, t=1,3)
                    IF (str=="B") THEN 
                        posB = 1
                        PRINT *, "La supernova ", "SN"//names(n), " ha la seconda colonna come colonna della magnitudine 'B'"
                    END IF
                ELSE
                    READ(11,*) tempMJD, (tempB, t=1,posB), tempBerr
                    IF (tempB.lt.99.d0) THEN !Salva i valori letti solamente se la magnitudine 'B' è minore di 99.d0
                        sn_mat(j,k) = tempMJD
                        sn_mat(j,k+1) = tempB
                        sn_mat(j,k+2) = tempBerr
                        j = j + 1
                    END IF
                END IF                
            END DO
            CLOSE(11)
            k = k + 3
        END DO
        PRINT *, "Salvataggio dei dati completato con successo."
        PRINT *, ""
    END SUBROUTINE sn_buildmat

    SUBROUTINE mask_buildmat(n_sn, dim1, dim2, sn_mat, mask, valid_sn) !Crea un array logico 'mask' contenente quali supernovae sono monotone e quali no
        INTEGER, INTENT(IN) :: dim1, dim2, n_sn
        REAL*8, DIMENSION(dim1, dim2), INTENT(IN) :: sn_mat
        INTEGER, INTENT(OUT) :: valid_sn
        LOGICAL, ALLOCATABLE, INTENT(OUT) :: mask(:)

        INTEGER :: i, j

        ALLOCATE(mask(n_sn))

        j=1
        valid_sn = 0
        PRINT *, "Controllo della monotonia delle curve di luce in corso..."
        DO i=2,dim2,3
            mask(j) = monotone_check(sn_mat(:,i),dim1)
            IF(mask(j)) valid_sn = valid_sn + 1 !Conta il numero di curve di luce non monotone
            j = j + 1
        END DO
        PRINT *, "Controllo della monotonia delle curve di luce completato con successo."
    END SUBROUTINE mask_buildmat

    SUBROUTINE export_mat(mat, names, size1_mat, size2_mat, size_names, filename, opt) !Si occupa della scrittura dei vari file di output
        INTEGER, INTENT(IN) :: size1_mat, size2_mat, size_names, opt
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(LEN=8), DIMENSION(size_names), INTENT(IN) :: names
        REAL*8, INTENT(IN) :: mat(size1_mat,size2_mat)

        INTEGER :: i, j
        CHARACTER(LEN=4) :: extension
        CHARACTER(LEN=10), DIMENSION(24) :: labels
        CHARACTER(LEN=8), DIMENSION(size_names*11) :: new_names
        CHARACTER(LEN=10), DIMENSION(size_names*11) :: new_labels

        extension = '.dat'
        new_names = ''
        new_labels = ''
        labels = ['names     ', 'zCMB      ', 'JD_max    ', '+/-       ', 'b_max     ', '+/-       ', 'b15       ', &
                '+/-       ', 'Delta_m15 ', '+/-       ', 'B_max     ', '+/-       ', 'mu        ', '+/-       ', &
                'dL        ', '+/-       ', 'H0        ', '+/-       ', 'H0_ma     ', '+/-       ', 'H0_mp     ', &
                '+/-       ', 'H0_mc     ', '+/-       ']

        j = 1
        OPEN(11, file = filename//extension)
        SELECT CASE (opt)
            CASE(0) !Opzione che scrive tutti i calcoli fatti in un unico file
                PRINT *, "Scrittura dei dati di tutte le supernovae nel file '"//filename//extension//"' in corso..."
                WRITE(11,'(1x, *(g0, ", "))') (labels(j), j=1,SIZE(labels))
                DO i = 1,SIZE(mat,1)
                    WRITE(11,'(1x, *(g0, ", "))') names(i), (mat(i,j), j=1,SIZE(mat,2))
                END DO

            CASE(1) !Opzione che scrive i calcoli fatti della spline per le supernovae selezionaate
                PRINT *, "Scrittura dei dati della spline nel file '"//filename//extension//&
                        "' per le supernovae richieste in corso..."
                DO i = 1, SIZE(names)
                    new_names(j) = names(i)
                    new_names(j+1) = ''
                    new_labels(j) = labels(3)
                    new_labels(j+1) = labels(5)
                    j = j + 2
                END DO
    
                WRITE(11,'(1x, *(g0, ", "))') "#"//new_names(1), (new_names(j), j=2,(SIZE(names)*2-1))
                WRITE(11,'(1x, *(g0, ", "))') "#"//new_labels(1), (new_labels(j), j=2,(SIZE(names)*2))
                DO i = 1,SIZE(mat,1)
                    WRITE(11,'(1x, *(g0, ", "))') (mat(i,j), j=1,SIZE(mat,2))
                END DO

            CASE(2) !Opzione che scrive i dati generati casualmente per le supernovae selezionate
                PRINT *, "Scrittura dei dati delle distribuzioni nel file '"//filename//extension//&
                        "' per le supernovae richieste in corso..."
                DO i = 1, SIZE(names)
                    new_names(j) = names(i)
                    new_names(j+1:j+10) = ''
                    new_labels(j:j+10) = labels(3:SIZE(labels):2)
                    j = j + 11
                END DO
    
                WRITE(11,'(1x, *(g0, ", "))') "#"//new_names(1), (new_names(j), j=2,SIZE(new_names))
                WRITE(11,'(1x, *(g0, ", "))') "#"//new_labels(1), (new_labels(j), j=2,SIZE(new_labels))
                DO i = 1,SIZE(mat,1)
                    WRITE(11,'(1x, *(g0, ", "))') (mat(i,j), j=1,SIZE(mat,2))
                END DO

        END SELECT
        CLOSE(11)
        PRINT *, "Scrittura del file completata con successo."
    END SUBROUTINE export_mat

END MODULE save_file

MODULE spline_interp
!Insieme di funzioni e subroutine necessarie per l'interpolazione tramite spline cubica
IMPLICIT NONE

    CONTAINS
    
    REAL*8 FUNCTION fun_e(arr,pos,n)
        INTEGER :: pos, n
        REAL*8 :: arr(0:n)

        fun_e = arr(pos) - arr(pos-1)
    END FUNCTION fun_e

    REAL*8 FUNCTION fun_r(arr,pos,n)
        INTEGER :: pos, n
        REAL*8 :: arr(0:n)

        fun_r = 2.*(arr(pos+1) - arr(pos-1))
    END FUNCTION fun_r

    REAL*8 FUNCTION fun_g(arr,pos,n)
        INTEGER :: pos, n
        REAL*8 :: arr(0:n)

        fun_g = arr(pos+1) - arr(pos)
    END FUNCTION fun_g

    REAL*8 FUNCTION fun_c(arr,fun,pos,n)
        INTEGER :: pos, n
        REAL*8 :: arr(0:n), fun(0:n)

        fun_c = (6/(arr(pos+1) - arr(pos)))*(fun(pos+1) - fun(pos)) + (6/(arr(pos) - arr(pos-1)))*(fun(pos-1) - fun(pos))
    END FUNCTION fun_c

    SUBROUTINE spline_buildmat(a,c,x,f,n) !Costruisce la matrice e i vettori necessari per il calcolo della spline
        INTEGER, INTENT(IN) :: n

        INTEGER :: i
        REAL*8 :: a(0:n,0:n), c(0:n), x(0:n), f(0:n)
        !'a' matrice tri-diagonale (n+1,n+1), 'c' vettore (n+1) dei termini noti, 'f' vettore (n+1) delle magnitudini osservate, 'x' vettore (n+1) dei giorni delle osservazioni
    
        a = 0.d0
        c = 0.d0
        DO i = 0,n
            IF(i == 0) THEN
                a(0,0) = 1.d0
            ELSE IF(i == n) THEN
                a(n,n) = 1.d0
            ELSE
                a(i,i-1) = fun_e(x,i,n)
                a(i,i) = fun_r(x,i,n)
                a(i,i+1) = fun_g(x,i,n)
    
                c(i) = fun_c(x,f,i,n)
            END IF
        END DO
    END SUBROUTINE spline_buildmat
    
    SUBROUTINE gauss(a,c,f2,n) !Risolve con Gauss il sistema richiesto dalla spline
        INTEGER :: n, i, j, k
        REAL*8, DIMENSION(n) :: c, f2 !'c' è il vettore dei termini noti e 'f2' è il vettore delle soluzioni
        REAL*8, DIMENSION(n,n) :: a !'a' è la matrice del sistema
        REAL*8 :: fakt, summa
        
        DO i = 1,n-1 !Sceglie la variabile da eliminare
            DO j = i+1, n !Sceglie la riga su cui eliminare
                fakt = a(j,i)/a(i,i)
                DO k = 1,n
                    a(j,k) = a(j,k) - a(i,k)*fakt
                END DO
                c(j) = c(j) - c(i)*fakt
            END DO
        END DO
    
        f2(n) = c(n)/a(n,n)
        DO i = n-1,1,-1
            summa = 0.d0
            DO j = i+1,n
                summa = summa + a(i,j)*f2(j)
            END DO
            f2(i) = (c(i) - summa)/a(i,i)
        END DO
    END SUBROUTINE gauss

    SUBROUTINE peak_finder(npoints, ndata, x, f, f2, pos_max, t_max, m_max, m_15, spline_mat, opt) !Trova il massimo della curva di luce tramite spline
        !'x' è il vettore dei giorni JD, 'f' è il vettore delle magnitudini osservate e 'f2' il vettore delle soluzioni del sistema
        INTEGER, INTENT(IN) :: ndata, pos_max, opt
        REAL*8, INTENT(IN) :: npoints
        REAL*8, DIMENSION(0:ndata-1), INTENT(IN) :: x, f, f2
        REAL*8, DIMENSION(INT(npoints),2), INTENT(OUT) :: spline_mat
        REAL*8, INTENT(OUT) :: t_max, m_max, m_15

        INTEGER :: i, j, max_cicle
        REAL*8 :: xin, fout, interval

        j = 1
        fout = 99.0d0
        m_max = 99.0d0

        !Calcola i punti in cui calcolare la spline
        IF (opt == 1) THEN 
            interval = (x(ndata-1) - x(0))/npoints !Se si tratta di una supernova richiesta, assume come intervallo l'intera curva di luce e ne calcola il passo
            xin = x(0)+interval
            max_cicle = ndata-1
        ELSE 
            interval = (x(pos_max) - x(pos_max-2))/npoints !Se non si tratta di una supernova richiesta, assume come intervallo solamente l'intorno del massimo e ne calcola il massimo
            xin = x(pos_max-2)+interval
            max_cicle = pos_max
        END IF
        
        DO WHILE (xin.lt.x(max_cicle))
            IF (opt == 1) THEN !Avanza seguendo il passo scelto sopra, partendo dal punto adeguato
                xin = x(0)+interval*j
            ELSE
                xin = x(pos_max-2)+interval*j
            END IF

            DO i = 1, ndata-1
                IF (xin.gt.x(i-1) .and. xin.le.x(i)) THEN !Trova la posizione del punto selezionato e applica la formula della spline
                    fout = (f2(i-1) * ((x(i) - xin)**3) + f2(i) * ((xin - x(i-1))**3)) / & 
                        (6.0d0*(x(i)-x(i-1))) + &
                        (f(i-1) / (x(i)-x(i-1)) - &
                        f2(i-1)*(x(i)-x(i-1)) / 6.0d0)*(x(i)-xin) + &
                        (f(i) / (x(i)-x(i-1)) - &
                        f2(i)*(x(i)-x(i-1)) / 6.0d0)*(xin - x(i-1))
                END IF
            END DO

            IF (opt == 1) THEN !Per le supernovae richieste vengono salvata tutti i valori del tempo e della magnitudine calcolati dalla spline
                spline_mat(j,1) = xin
                spline_mat(j,2) = fout
            END IF

            IF (fout.lt.m_max) THEN !Trova la magnitudine al massimo ed il relativo tempo
                m_max = fout
                t_max = xin
            END IF
            j = j + 1
        END DO

        xin = t_max + 15
        DO i = MAX(pos_max-2,1), ndata-1
            IF (xin.gt.x(i-1) .and. xin.lt.x(i)) THEN !Applica la spline a 15 giorni dopo il massimo
                m_15 = (f2(i-1) * ((x(i) - xin)**3) + f2(i) * ((xin - x(i-1))**3)) / & 
                    (6.0d0*(x(i)-x(i-1))) + &
                    (f(i-1) / (x(i)-x(i-1)) - &
                    f2(i-1)*(x(i)-x(i-1)) / 6.0d0)*(x(i)-xin) + &
                    (f(i) / (x(i)-x(i-1)) - &
                    f2(i)*(x(i)-x(i-1)) / 6.0d0)*(xin - x(i-1))
            END IF
        END DO
    END SUBROUTINE peak_finder

    SUBROUTINE spline(sn_data, ndati, pos_lastval, JD_max, b_max, b_15, spline_mat, npoints, opt) !Gestisce le varie subroutine che si occupano del calcolo della spline
        INTEGER, INTENT(IN) :: ndati, npoints, opt
        REAL*8, INTENT(IN) :: sn_data(ndati,2)
        INTEGER, INTENT(INOUT) :: pos_lastval
        REAL*8, DIMENSION(npoints,2), INTENT(OUT) :: spline_mat
        REAL*8, INTENT(OUT) :: JD_max, b_max, b_15

        INTEGER :: pos_max
        REAL*8, ALLOCATABLE :: a(:,:), f2(:), c(:)

        IF (pos_lastval == 0) THEN !Trova la posizione dell'ultimo valore
            DO pos_lastval = ndati, 1, -1
                IF (sn_data(pos_lastval,2).lt.99.0d0) EXIT
            END DO
        END IF

        DO pos_max = 1, pos_lastval !Trova la posizione della magnitudine più piccola dei dati osservati
            IF (sn_data(pos_max,2) == MINVAL(sn_data(1:pos_lastval,2))) EXIT
        END DO


        ALLOCATE(a(pos_lastval,pos_lastval),f2(pos_lastval),c(pos_lastval))

        f2=0.d0
        CALL spline_buildmat(a,c,sn_data(1:pos_lastval,1),sn_data(1:pos_lastval,2),pos_lastval-1)
        CALL gauss(a,c,f2,pos_lastval)
        CALL peak_finder(REAL(npoints,8),pos_lastval,sn_data(1:pos_lastval,1),sn_data(1:pos_lastval,2),&
            f2,pos_max,JD_max,b_max,b_15,spline_mat,opt)
    END SUBROUTINE spline

END MODULE spline_interp

MODULE quantities
!Insieme di funzioni che calcolano tutte le quantità richieste (tranne le costanti di Hubble) a partire dai dati ottenuti dalla spline
IMPLICIT NONE

    CONTAINS

    REAL*8 FUNCTION delta_m15(b_max, b15) !Calcola la differenza di magnitudine a 15 giorni dal massimo
        REAL*8, INTENT(IN) :: b_max, b15

        delta_m15 = b15 - b_max
    END FUNCTION delta_m15

    REAL*8 FUNCTION correction(b_obs, BV, Rv) !Corregge per l'estinzione la magnitudine al massimo
        REAL*8, INTENT(IN) :: b_obs, BV, Rv

        correction = b_obs - (Rv + 1.d0)*BV
    END FUNCTION correction

    REAL*8 FUNCTION M_abs(m15) !Calcola la magnitudine assoluta
        REAL*8, INTENT(IN) :: m15

        M_abs = -19.258 + .784*(m15 - 1.1)
    END FUNCTION M_abs

    REAL*8 FUNCTION distance_module(b_obs,B_abs) !Calcola il modulo di distanza
        REAL*8, INTENT(IN) :: b_obs, B_abs

        distance_module = b_obs - B_abs
    END FUNCTION distance_module

    REAL*8 FUNCTION luminosity_distance(mu) !Calcola la distanza di luminosità
        REAL*8, INTENT(IN) :: mu

        luminosity_distance = 10**((mu + 5.d0)*.2d0 - 6.d0)
    END FUNCTION luminosity_distance
    
END MODULE quantities

MODULE hubble
!Insieme di funzioni e subroutine necessarie a calcolare le costanti di Hubble per i vari modelli considerati
IMPLICIT NONE

    CONTAINS

    REAL*8 FUNCTION integrand(z, Omega_M, Omega_Lambda) !Funzione integranda di cui si vuole calcolare l'integrale
        REAL*8 :: z, Omega_M, Omega_Lambda

        integrand = 1/SQRT(Omega_M*(1 + z)**3 + Omega_Lambda + (1 - Omega_M - Omega_Lambda)*(1 + z)**2)
    END FUNCTION integrand

    SUBROUTINE simpson1_3(a, b, n, Omega_M, Omega_Lambda, res) !Implementazione del metodo di risoluzione numerica dell'integrale
        INTEGER, INTENT(IN) :: n
        REAL*8, INTENT(IN) :: a, b, Omega_M, Omega_Lambda
        REAL*8, INTENT(OUT) :: res

        INTEGER :: i
        REAL*8 :: h, sum1, sum2

        h = (b-a)/n !Lunghezza degli n sottointervalli
        sum1 = 0.0d0
        sum2 = 0.0d0

        DO i = 1, n - 1, 2 !Calcola la prima sommatoria di Simpson 1/3
            sum1 = sum1 + integrand(a + i*h, Omega_M, Omega_Lambda)
        END DO

        DO i = 2, n - 2, 2 !Calcola la seconda sommatoria di Simpson 1/3
            sum2 = sum2 + integrand(a + i*h, Omega_M, Omega_Lambda)
        END DO

        res = h/3.0d0*(integrand(a, Omega_M, Omega_Lambda) + 4*sum1 + 2*sum2 + integrand(b, Omega_M, Omega_Lambda)) !Trova il valore numerico dell'integrale richiesto

    END SUBROUTINE simpson1_3

    SUBROUTINE hubble_costant_int(zCMB, dL, H0, H0_ma, H0_mp, H0_mc) !Risolve gli integrali delle costanti di Hubble delle singole supernovae per i vari modelli
        REAL*8, INTENT(IN) :: zCMB, dL
        REAL*8, INTENT(OUT) :: H0, H0_ma, H0_mp, H0_mc

        INTEGER, PARAMETER :: n = 1000 !Numero di sottointervalli di cui dividere l'intervallo in Simpson 1/3
        REAL*8, PARAMETER :: Omega_M = 0.3d0, Omega_Lambda = 0.7d0
        REAL*8, PARAMETER :: Omega_M_ma = 0.3d0, Omega_Lambda_ma = 0.d0
        REAL*8, PARAMETER :: Omega_M_mp = 1.d0, Omega_Lambda_mp = 0.d0
        REAL*8, PARAMETER :: Omega_M_mc = 2.d0, Omega_Lambda_mc = 0.d0
        REAL*8, PARAMETER :: c_light = 299792.458d0
        REAL*8 :: a, b, res1, res2, res3, res4

        a = 0.d0
        b = zCMB

        CALL simpson1_3(a, b, n, Omega_M, Omega_Lambda, res1)
        H0 = c_light*(1+zCMB)*res1/dL !Costante di Hubble del modello di riferimento
        CALL simpson1_3(a, b, n, Omega_M_ma, Omega_Lambda_ma, res2)
        H0_ma = c_light*(1+zCMB)*res2/dL !Costante di Hubble del modello aperto
        CALL simpson1_3(a, b, n, Omega_M_mp, Omega_Lambda_mp, res3)
        H0_mp = c_light*(1+zCMB)*res3/dL !Costante di Hubble del modello piatto
        CALL simpson1_3(a, b, n, Omega_M_mc, Omega_Lambda_mc, res4)
        H0_mc = c_light*(1+zCMB)*res4/dL !Costante di Hubble del modello chiuso

    END SUBROUTINE hubble_costant_int

    SUBROUTINE hubble_costant(H0_arr, err_arr, dim, H0_est, H0_err_est) !Calcola la media pesata della costante di Hubble
        INTEGER, INTENT(IN) :: dim
        REAL*8, DIMENSION(dim), INTENT(IN) :: H0_arr, err_arr
        REAL*8, INTENT(OUT) :: H0_est, H0_err_est

        INTEGER :: i
        REAL*8 :: num, den

        num = 0.0d0
        den = 0.0d0
        DO i = 1, dim
            num = num + (H0_arr(i)/err_arr(i)**2)
            den = den + 1/(err_arr(i)**(2))
        END DO

        H0_est = num/den !Media pesata della costante di Hubble
        H0_err_est = SQRT(1/den) !Errore sulla media pesata della costante di Hubble
    

    END SUBROUTINE hubble_costant

END MODULE hubble

MODULE errors
!Insieme di funzioni e subroutine necessarie per il calcolo degli errori
IMPLICIT NONE

    INTEGER, PRIVATE :: jran = 2

    CONTAINS

    REAL*8 FUNCTION lingen() !Genera uniformemente numeri che vanno da 0 ad 1
        INTEGER, PARAMETER :: mran = 259200
        INTEGER, PARAMETER :: cran = 54773
        INTEGER, PARAMETER :: aran = 7141

        jran = MOD(aran*jran + cran, mran)
        lingen = FLOAT(jran)/FLOAT(mran)

    END FUNCTION lingen

    REAL*8 FUNCTION percentile(arr,dim) !Calcola gli errori tramite 16-esimo e 84-esimo percentile
        INTEGER, INTENT(IN) :: dim
        REAL*8, INTENT(INOUT) :: arr(dim)

        CALL sorting(arr,dim)

        percentile = arr(NINT(dim*0.84d0)) - arr(NINT(dim*0.16d0))

    END FUNCTION percentile

    SUBROUTINE set_seed(j) !Imposta il seme della generazione casuale dei numeri
        INTEGER, INTENT(IN) :: j

        jran = j
    END SUBROUTINE set_seed

    SUBROUTINE box_muller(mat_data,dim,pos_lastval,new_mat) !Genera numeri con una distribuzione gaussiana
        INTEGER, INTENT(IN) :: dim, pos_lastval
        REAL*8, DIMENSION(dim,3), INTENT(IN) :: mat_data
        REAL*8, DIMENSION(dim,2), INTENT(OUT) :: new_mat

        INTEGER :: i
        REAL*8, PARAMETER :: pi = ACOS(-1.0d0)
        REAL*8 :: x1, x2, y1, y2
        
        new_mat = 99.0d0
        new_mat(:,1) = mat_data(:,1)
        i = 1
        DO i=1,pos_lastval,2
            x1 = lingen()
            x2 = lingen()
    
            !Numeri con una distribuzione gaussiana di media nulla e varianza unitaria
            y1 = SQRT(-2.0d0*LOG(x1))*COS(2*pi*x2)
            y2 = SQRT(-2.0d0*LOG(x1))*SIN(2*pi*x2)
    
            !Numeri con una distribuzione gaussiana di media il valore osservato e varianza l'errore del valore osservato
            new_mat(i,2) = y1*mat_data(i,3) + mat_data(i,2)
            IF(MOD(pos_lastval,2) == 0 .or. i /= pos_lastval) new_mat(i+1,2) = y2*mat_data(i+1,3) + mat_data(i+1,2)
    
        END DO
    
    END SUBROUTINE box_muller

    SUBROUTINE sorting(arr,dim) !Ordina l'array dato in ordine crescente
        INTEGER, INTENT(IN) :: dim
        REAL*8, INTENT(INOUT) :: arr(dim)

        INTEGER :: i, j, pos_min
        REAL*8 :: min_val
    

        DO i = 1, dim-1
            min_val = arr(i)
            pos_min = i
            DO j = i+1, dim
                IF (arr(j) .le. min_val) THEN
                    min_val = arr(j)
                    pos_min = j
                END IF
            END DO

            arr(pos_min)=arr(i)
            arr(i)=min_val

        END DO
    END SUBROUTINE sorting
    
END MODULE errors

PROGRAM main
USE spline_interp
USE save_file
USE quantities
USE hubble
USE errors
IMPLICIT NONE
    INTEGER, PARAMETER :: nran = 1000 !Numero di curve di luce generate casualmente per ogni supernova
    INTEGER, ALLOCATABLE, DIMENSION(:) :: valid_data
    INTEGER :: i, j, k, h1, h2, e1, e2, cicle !Variabili di ciclo
    INTEGER :: n_sn, valid_sn, seed, npoints, npoints_max, opt

    REAL*8, ALLOCATABLE, DIMENSION(:,:) :: newtable, sn_mat, results_mat, temp_mat, err_mat, temp, spline_mat, distribution_mat
    REAL*8 :: H0, err_H0, H0_ma, err_H0_ma, H0_mp, err_H0_mp, H0_mc, err_H0_mc !valori delle costanti di hubble ed i relativi errori

    CHARACTER(LEN=*), PARAMETER :: newtable_name = 'newtable2.txt'
    CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: valid_names, valid_export
    CHARACTER(LEN=6), ALLOCATABLE, DIMENSION(:) :: names, export_sn

    LOGICAL, PARAMETER :: view = .false. !Se 'view = .true.' verranno visualizzati a schermo i nomi delle supernovae le cui curve di luce sono monotone
    LOGICAL, PARAMETER :: skip = .false. !Se 'skip = .true.' non verranno estratti i dati come la curva di luce e la distribuzione delle quantità per le supernovae richieste
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: mask 

    !Creazione degli array dei dati dei file
    CALL newtable_buildmat(newtable_name, n_sn, newtable, names)
    CALL sn_buildmat(n_sn, names, sn_mat)
    CALL mask_buildmat(n_sn, SIZE(sn_mat,1), SIZE(sn_mat,2), sn_mat, mask, valid_sn)

    PRINT *, "Sulle", n_sn,"supernovae lette,", valid_sn,"supernovae sono non monotone."
    IF (view) THEN
        PRINT *, "Le supernovae monotone sono le seguenti:"
        DO i = 1, n_sn
            IF (mask(i) .eqv. .false.) PRINT *, 'SN'//names(i)
        END DO
    ELSE
        PRINT *, "Per vedere quali sono le supernovae monotone impostare 'view = .true.' all'interno del codice."
    END IF
    PRINT *, ""

    !Alloca gli array dove verranno salvati i risultati dei calcoli
    ALLOCATE(valid_data(n_sn), results_mat(valid_sn,23), valid_names(valid_sn))
    ALLOCATE(temp_mat(SIZE(sn_mat,1),2), err_mat(nran,11))
    valid_data = 0
    results_mat = 99.0d0

    !Vengono definite le supernovae per cui è richiesta l'estrazione dei dati e si verifica che tra queste non ci siano curve di luce monotone
    ALLOCATE(export_sn(5))
    export_sn = ['2004dt', '2005A ', '2004gc', '2005el', '2006gt'] !Il numero di elementi deve coincidere con il numero con cui è stata allocata
    IF (skip .eqv. .false.) PRINT *, "Sono stati richiesti i dati della spline e della distribuzione delle quantita' per",&
                                    SIZE(export_sn), "supernovae."

    h1 = 0
    j = 0
    IF (skip .eqv. .false.) THEN
        DO i = 1, n_sn
            IF (ANY(export_sn == names(i))) THEN
                IF (mask(i) .eqv. .false.) THEN 
                    PRINT *, "La supernova ", 'SN'//names(i), &
                        " di cui sono stati richiesti i dati si tratta di una supernova monotona."
                    h1 = h1 + 1
                END IF
                j = j + 1
            END IF
        END DO
        IF(h1 == 1) PRINT *, "Per tale supernova non e' dunque possibile estrarre i dati richiesti."
        IF(h1 .gt. 1) PRINT *, "Per tali supernovae non e' dunque possibile estrarre i dati richiesti."
    END IF

    npoints_max = 2000 !Numero di punti in cui calcolare la spline per le supernovae di cui bisogna fare il grafico della curva di luce
    ALLOCATE(temp(npoints_max, 2), spline_mat(1,1), distribution_mat(1,1), valid_export(1))
    IF (j == 0) THEN
        npoints_max = -999 !Assume questo valore se non deve calcolare la spline sulle supernovae selezionate
        IF (skip) THEN
            PRINT *, "Non verranno estratti i dati della spline e delle distribuzioni per le supernovae richieste."
            PRINT *, "Impostare 'skip = .false.' all'interno del codice per estrarli."
        ELSE
            PRINT *, "Non sono state trovate le supernovae di cui sono stati richiesti i dati."
        END IF
    ELSE
        IF (h1 .lt. SIZE(export_sn)) THEN !Se le supernovae selezionate non sono tutte monotone alloca gli array e li inizializza
            DEALLOCATE(spline_mat, distribution_mat, valid_export)
            ALLOCATE(valid_export(SIZE(export_sn)-h1))
            ALLOCATE(spline_mat(npoints_max, (SIZE(export_sn)-h1)*2), distribution_mat(nran,(SIZE(export_sn)-h1)*11))
            spline_mat = 99.0d0
            distribution_mat = 99.0d0
        ELSE
            npoints_max = -999
            PRINT *, "Le supernovae di cui sono stati richiesti i dati sono tutte monotone."
        END IF
    END IF
    PRINT *, ""

    k = 1
    j = 1
    h1 = 1
    h2 = 1
    seed = 1
    PRINT *, "Calcolo dei dati in corso..."
    DO i = 1, n_sn
        opt = 0
        npoints = 1000

        IF (mask(i)) THEN !Per le supernovae richieste calcola la spline su un numero di punti superiore
            IF ((skip .eqv. .false.) .and. (ANY(export_sn == names(i)))) THEN
                valid_export(NINT(h1/2.d0)) = 'SN'//names(i)
                npoints = npoints_max
                opt = 1
            END IF
            valid_names(k) = 'SN'//names(i)

            PRINT *, "Inizio dei calcoli della supernova SN"//names(i)//" in posizione", k, "su", valid_sn
            !Salva i calcoli con i dati letti dai file nelle corrispettive posizioni della matrice
            CALL spline(sn_mat(:,[j,j+1]),SIZE(sn_mat,1),valid_data(i),&
                        results_mat(k,2),results_mat(k,4),results_mat(k,6),temp,npoints,opt)
            results_mat(k,1) = newtable(i,1)
            results_mat(k,8) = delta_m15(results_mat(k,4), results_mat(k,6))
            results_mat(k,10) = M_abs(results_mat(k,8))
            results_mat(k,12) = distance_module(correction(results_mat(k,4),newtable(i,2),newtable(i,3)),results_mat(k,10))
            results_mat(k,14) = luminosity_distance(results_mat(k,12))
            CALL hubble_costant_int(results_mat(k,1), results_mat(k,14), results_mat(k,16), results_mat(k,18), &
                results_mat(k,20), results_mat(k,22))

            IF (opt == 1) THEN !Salva nella matrice 'spline_mat' i dati della spline delle supernovae richieste
                spline_mat(:,h1) = temp(:,1)
                spline_mat(:,h1+1) = temp(:,2)
                h1 = h1 + 2
                opt = 2
            END IF

            cicle = 0
            err_mat = 99.0d0
            DO WHILE (cicle .lt. nran)
                CALL set_seed(seed)
                CALL box_muller(sn_mat(:,j:j+2),SIZE(sn_mat,1),valid_data(i),temp_mat) !Genera casualmente le 'nran' curve di luce richieste per il calcolo degli errori
                IF (monotone_check(temp_mat(:,2),SIZE(temp_mat,1))) THEN !Calcola le quantità richieste solamente se la curva di luce generata casualmente non è monotona
                    cicle = cicle + 1
                    CALL spline(temp_mat,SIZE(temp_mat,1),valid_data(i),&
                                err_mat(cicle,1),err_mat(cicle,2),err_mat(cicle,3),temp,npoints,opt)
                    err_mat(cicle,4) = delta_m15(err_mat(cicle,2), err_mat(cicle,3))
                    err_mat(cicle,5) = M_abs(err_mat(cicle,4))
                    err_mat(cicle,6) = distance_module(correction(err_mat(cicle,2),newtable(i,2),newtable(i,3)),err_mat(cicle,5))
                    err_mat(cicle,7) = luminosity_distance(err_mat(cicle,6))
                    CALL hubble_costant_int(results_mat(k,1), err_mat(cicle,7), err_mat(cicle,8), err_mat(cicle,9), &
                        err_mat(cicle,10), err_mat(cicle,11))

                END IF
                seed = seed + 1 !Cambia ad ogni ciclo il seme della generazione
            END DO

            IF (opt == 2) THEN !Salva nella matrice 'distribution_mat' i dati della distribuzione delle supernovae richieste
                distribution_mat(:,h2:h2+10) = err_mat(:,:)
                h2 = h2 + 11
            END IF

            e2 = 1
            DO e1 = 3, 23, 2
                results_mat(k,e1) = percentile(err_mat(:,e2),SIZE(err_mat,1)) !Tramite i dati ottenuti sopra calcola l'errore da associare alle quantità calcolate
                e2 = e2 + 1
            END DO
            k = k + 1
        END IF
        j = j + 3
    END DO
    PRINT *, "Tutti i calcoli sono stati completati con successo."
    PRINT *, ""

    IF (npoints_max .gt. 0) THEN !Se sono state richieste delle supernovae la cui curva di luce non è monotona, allora scrive i dati ottenuti su file
        CALL export_mat(spline_mat,valid_export,&
                        SIZE(spline_mat,1),SIZE(spline_mat,2), SIZE(valid_export),"spline_data",1)
        CALL export_mat(distribution_mat,valid_export,&
                        SIZE(distribution_mat,1),SIZE(distribution_mat,2), SIZE(valid_export),"distribution_data",2)
    END IF

    CALL export_mat(results_mat,valid_names,SIZE(results_mat,1),SIZE(results_mat,2), SIZE(valid_names),"data",0) !Salva i dati di tutte le supernovae in un unico file

    !Calcola la media pesata delle costanti di Hubble per i vari modelli considerati
    PRINT *, ""
    PRINT *, "Calcolo delle costanti di Hubble per i vari modelli in corso..."
    CALL hubble_costant(results_mat(:,16),results_mat(:,17),SIZE(results_mat,1),H0,err_H0)
    CALL hubble_costant(results_mat(:,18),results_mat(:,19),SIZE(results_mat,1),H0_ma,err_H0_ma)
    CALL hubble_costant(results_mat(:,20),results_mat(:,21),SIZE(results_mat,1),H0_mp,err_H0_mp)
    CALL hubble_costant(results_mat(:,22),results_mat(:,23),SIZE(results_mat,1),H0_mc,err_H0_mc)

    WRITE(*,*) "Il valore stimato della costante di Hubble con il modello di riferimento e'", H0, "con errore", err_H0
    WRITE(*,*) "Il valore stimato della costante di Hubble con il modello aperto e'        ", H0_ma, "con errore", err_H0_ma
    WRITE(*,*) "Il valore stimato della costante di Hubble con il modello piatto e'        ", H0_mp, "con errore", err_H0_mp
    WRITE(*,*) "Il valore stimato della costante di Hubble con il modello chiuso e'        ", H0_mc, "con errore", err_H0_mc
    
END PROGRAM main

!!! STRUTTURA DELLA MATRICE results_mat !!!
!   results_mat(:,1)        = redshift                                                      |zCMB
!   results_mat(:,[2,3])    = tempo giuliano del massimo e relativo errore                  |JD_max
!   results_mat(:,[4,5])    = magnitudine 'b' al massimo e relativo errore                  |b_max
!   results_mat(:,[6,7])    = magnitudine 'b' dopo 15 giorni e relativo errore              |b15
!   results_mat(:,[8,9])    = differenza tra 'b15' e 'b_max' e relativo errore              |Delta_m15
!   results_mat(:,[10,11])  = magnitudine 'B' al massimo e relativo errore                  |B_max
!   results_mat(:,[12,13])  = modulo di distanza e relativo errore                          |mu
!   results_mat(:,[14,15])  = distanza di luminosità e relativo errore                      |dL
!   results_mat(:,[16,17])  = costante di hubble e relativo errore                          |H0
!   results_mat(:,[18,19])  = costante di hubble per il modello aperto e relativo errore    |H0_ma
!   results_mat(:,[20,21])  = costante di hubble per il modello piatto e relativo errore    |H0_mp
!   results_mat(:,[22,23])  = costante di hubble per il modello chiuso e relativo errore    |H0_mc