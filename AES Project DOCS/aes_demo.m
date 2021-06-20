%AES_DEMO  Demonstration of AES-components.
%
%   AES_DEMO
%   runs a demonstration of all components of 
%   the Advanced Encryption Standard (AES) toolbox.
%
%   In the initialization step the S-boxes, the round constants,
%   and the polynomial matrices are created and
%   an example cipher key is expanded into 
%   the round key schedule.
%   Step two and three finally convert 
%   an example plaintext to ciphertext and back to plaintext.

%   Copyright 2001-2005, J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de

%   Version 1.0     30.05.2001

% Initialization
clc
disp('START')
ciphertext1=zeros(1,16);
[s_box, inv_s_box, w, poly_mat, inv_poly_mat] = aes_init;

% RAW AES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=imread('girl.jpg');
% [row,col,dim]=size(A);
% B=A;
% %B=rgb2gray(A);
% C=dec2hex(B);
% disp(C)
% D=zeros(row,col);
% E=zeros(row,col);
% [r, c] =size(C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0:2570;
left=mod(a.^2,3049);
right=mod(a.^3-a+188,3049);
points=[];
for i=1:length(right)
    I=find(left==right(i));
    for j=1:length(I)
        points=[points;a(i),a(I(j))];
    end
end
k=10;
store=1;
store1=1;
mat=imread('camera_128_gray.jpg');
orig=mat
[r c]=size(mat);
q=zeros(1,r*c);
q1=zeros(1,r*c);
for l=1:r
    for m=1:c
        num=mat(l,m);
        value=10*int16(num);
        flag=0;
        for n=1:11
            if flag==1
                break;
            end
            compare=value+n-1;
            for o=1:2145
                if compare==points(o,1)
                    q(1,store)=fix(points(o,2)/255);
                    store=store+1;
                    mat(l,m)=mod(points(o,2),255);
                    flag=1;
                    q1(1,store1)=compare;
                    store1=store1+1;
                    break;
                end
            end
        end
    end
end
ecc_mat=mat
%mat
% imshow(mat)
% mat=imread('lena1.jpg');
% mat=rgb2gray(mat);
[m,n] = size(mat);
inputMat=mat
key = [34, 65, 98, 23, 56, 99, 14, 00, 89, 35, 12, 16, 67, 69, 25, 49];   %this is a 16 byte key
[row,col] = size(inputMat);
[m,n] = size(key);
for i=1:row
    for j=1:col
        mat(i,j) = s_box(inputMat(i,j)+1);   % s_box function returns the mapped s_box value 
    end
end
oneDimen = reshape(mat, [1,row*col]);  % reshape converts one matrix into another matrix of required dimensions 
x = mod(row*col,m*n);
if x~=0
    chunks = ((row*col)/(m*n))+1;
else
    chunks = (row*col)/(m*n);
end
for i=1:chunks
    encryShiftedKey = circshift(key, [0,i-1]);
    for j=1:m*n
        index = (i-1)*(m*n) + j;
        if index > row*col
            break;
        end
        oneDimen(index) = bitxor(oneDimen(index), encryShiftedKey(j));
    end
end
magicMatrix = magic(row);
for i=1:row
    for j=1:row
        magicMatrix(i,j) = oneDimen(magicMatrix(i,j));
    end
end
newmat=magicMatrix
[row col]=size(newmat);
B=uint8(newmat);
C=dec2hex(newmat);
disp(C)
D=zeros(row,col);
E=zeros(row,col);
[r, c] =size(C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define an arbitrary series of 16 plaintext bytes 
% in hexadecimal (string) representation
% The following two specific plaintexts are used as examples 
% in the AES-Specification (draft)
% plaintext_hex=array;
%
% ENCRYPTION
%
x=1;
y=1;
for i=1:16:r
%i=1;
    plaintext_hex = {strcat(C(i,1),C(i,2)) strcat(C(i+1,1),C(i+1,2)) strcat(C(i+2,1),C(i+2,2)) strcat(C(i+3,1),C(i+3,2)) strcat(C(i+4,1),C(i+4,2)) strcat(C(i+5,1),C(i+5,2)) strcat(C(i+6,1),C(i+6,2)) strcat(C(i+7,1),C(i+7,2)) ...
                     strcat(C(i+8,1),C(i+8,2)) strcat(C(i+9,1),C(i+9,2)) strcat(C(i+10,1),C(i+10,2)) strcat(C(i+11,1),C(i+11,2)) strcat(C(i+12,1),C(i+12,2)) strcat(C(i+13,1),C(i+13,2)) strcat(C(i+14,1),C(i+14,2)) strcat(C(i+15,1),C(i+15,2))};
    %plaintext_hex = {'32' '43' 'f6' 'a8' '88' '5a' '30' '8d' ...
    %                 '31' '31' '98' 'a2' 'e0' '37' '07' '34'};

    % Convert plaintext from hexadecimal (string) to decimal representation
    plaintext = hex2dec (plaintext_hex);

    % MAGIC MATRIX ENC
%     M=magic(4);
%     T=zeros(4,4);
%     for m=1:4
%         for n=1:4
%         temp = M(m,n);
%         T(m,n)=plaintext(temp);
%         end
%     end
%     T=T(:)';
%     plaintext=T;
    %    
    ciphertext = cipher (plaintext, w, s_box, poly_mat, 1);
    for z=1:16
        D(x,y)=ciphertext(1,z);
        x=x+1;
        if x>row
            x=1;
            y=y+1;
        end
        if y>col
            break;
        end
    end
end
    % Convert the ciphertext back to plaintext
    % using the expanded key, the inverse S-box, 
    % and the inverse polynomial transformation matrix
    disp(D)
%
%DECRYPTION
%
    Z=D;
    x=1;
    y=1;
    D=D(:)';
    [nr, nc]=size(D);
    for i=1:16:nc
        for j=1:16
            ciphertext1(1,j)=D(1,i+j-1);
        end
        re_plaintext = inv_cipher (ciphertext1, w, inv_s_box, inv_poly_mat, 1);
        % MAGIC MATRIX DECYP
        %conv replaintext to 4x4
%         T=zeros(1,16);
%         for m=1:4
%             for n=1:4
%         T(M(m,n))=replain(m,n);
%             end
%         end
%         re_plaintext=replain;    
        %
        for z=1:16
        E(x,y)=re_plaintext(1,z);
        x=x+1;
            if x>row
                x=1;
                y=y+1;
            end
            if y>col
                break;
            end
        end
    end 
    final=uint8(E);
    aes_op = uint8(E)
%     subplot(1,3,1),imshow(B);
%     title('original');
%     subplot(1,3,2),imshow(uint8(Z));
%     title('encrypted');
%     subplot(1,3,3),imshow(final);
%     title('decrypted');
%     disp('END')
 %%% MM DECRYP
magicMatrix=final
resultOneDMatrix = zeros(1,row*col);
dummyMagicMatrix = magic(row);
for i=1:row
    for j=1:row
        resultOneDMatrix(dummyMagicMatrix(i,j)) = magicMatrix(i,j);
    end
end
for i=1:chunks
    decryShiftedKey = circshift(key, [0,i-1]);
    for j=1:m*n
        index = (i-1)*(m*n) + j;
        if index > row*col
            break;
        end
        resultOneDMatrix(index) = bitxor(resultOneDMatrix(index), decryShiftedKey(j));
    end
end
result = reshape(resultOneDMatrix, [row,col]);
for i=1:row
    for j=1:col
        newMat(i,j) = inv_s_box(result(i,j)+1);   % inv_s_box function returns the inverse mapped s_box value 
    end
end

finalmat=newMat
mat=uint8(finalmat);
[r c]=size(mat);
%ecc decrypt
s=1;
z=1;
compare=0;
for i=1:r
    for j=1:c
           num2=int16(mat(i,j));
           num2=q(1,s)*255 + num2;
           s=s+1;
           compare= num2;
           for l=1:2145
               if compare == points(l,2) && q1(1,z)==points(l,1)
                     mat(i,j)=points(l,1)/10;
                     break;
               end    
           end
           z=z+1;
    end
end

final=mat
display(final);
imshow(final)