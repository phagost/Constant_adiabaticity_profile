function [Ix,Iy,Iz,Ip,Im] = Operators(N,I)
    %{
    This function is made for defining operators of N-spins system with the spin I
    input : the number of spins N; the spin value I
    result : Ix, Iy, Iz, Ip and Im operators
    %}

    %define an identity matrix
    E=eye(2*I+1);

    %define Zeeman basis for the one nucleus
    %I have to think how to make this array of arraies predefined
    b=cell(1,2*I+1);
    for k=1:(2*I+1)
        b{k}=diag(zeros(2*I+1));
        b{k}(k)=1;
    end

    %define operators for the one nucleus
    Ix0=zeros(2*I+1); Iy0=zeros(2*I+1);
    Iz0=zeros(2*I+1); Ip0=zeros(2*I+1);
    Im0=zeros(2*I+1);

    for k=1:(2*I+1)
        for j=1:(2*I+1)
            Iz0(k,j)=(b{k}'*b{j})*(-j+(I+1));

            if abs(j-1)<1
                Ip0(k,j)=0;
            else
                Ip0(k,j)=(b{k}'*b{j-1})*sqrt((I-(-j+(I+1)))*(I+(-j+(I+1))+1));
            end

            if abs(j+1)>(2*I+1)
                Im0(k,j)=0;
            else
                Im0(k,j)=(b{k}'*b{j+1})*sqrt((I+(-j+(I+1)))*(I-(-j+(I+1))+1));
            end
        end
    end

    Ix0=(Ip0+Im0)/2;
    Iy0=(Ip0-Im0)/(2*1i);

    %define operators for N nuclei
    Ix=cell(1,N); Iy=cell(1,N); Iz=cell(1,N);
    Ip=cell(1,N); Im=cell(1,N);
    for k=1:N
        for j=1:N
            if (k==1) && (j==1) 
                Ix{k}=Ix0;
                Iy{k}=Iy0;
                Iz{k}=Iz0;
                Ip{k}=Ip0;
                Im{k}=Im0;
            elseif (k~=1) && (j==1)
                Ix{k}=E;
                Iy{k}=E;
                Iz{k}=E;
                Ip{k}=E;
                Im{k}=E;
            elseif (k==j)
                Ix{k}=kron(Ix{k},Ix0);
                Iy{k}=kron(Iy{k},Iy0);
                Iz{k}=kron(Iz{k},Iz0);
                Ip{k}=kron(Ip{k},Ip0);
                Im{k}=kron(Im{k},Im0);
            else
                Ix{k}=kron(Ix{k},E);
                Iy{k}=kron(Iy{k},E);
                Iz{k}=kron(Iz{k},E);
                Ip{k}=kron(Ip{k},E);
                Im{k}=kron(Im{k},E);
            end
        end
    end
end

