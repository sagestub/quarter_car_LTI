
% ws_rankfunctions_kt_np.m  ver 1.0  by Tom Irvine

function[iwin,nrank,dm,drank]=...
    ws_rankfunctions_kt_np(rntrials,ym,vm,dm,em,im,cm,km,dskm,nspec,aunit,vunit,dunit,displacement_limit,iw,ew,dw,vw,aw,cw,kw,dskw)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


    yrank=zeros(rntrials,1);
	vrank=zeros(rntrials,1);
	drank=zeros(rntrials,1);
	erank=zeros(rntrials,1);
	irank=zeros(rntrials,1);
	crank=zeros(rntrials,1);    
    krank=zeros(rntrials,1);
    dskrank=zeros(rntrials,1);    
 
    pyrank=zeros(rntrials,1);
	pvrank=zeros(rntrials,1);
	pdrank=zeros(rntrials,1);
	perank=zeros(rntrials,1);
	pirank=zeros(rntrials,1);
    pcrank=zeros(rntrials,1);     
    pkrank=zeros(rntrials,1);    
    pdskrank=zeros(rntrials,1);   
    
	for i=1:rntrials
%
		yrank(i)=i;
		vrank(i)=i;
		drank(i)=i;
		erank(i)=i;
		irank(i)=i;
        crank(i)=i;        
        krank(i)=i;
        dskrank(i)=i;
%
    end
%
	for i=1:rntrials-1
%
		for j=i+1:rntrials
%
			if(ym(i)<ym(j))
				temp=ym(i);
                ym(i)=ym(j);  
				ym(j)=temp;
%
				itemp=yrank(i);
                yrank(i)=yrank(j);  
				yrank(j)=itemp;
            end    
			if(vm(i)<vm(j))
				temp=vm(i);
                vm(i)=vm(j);  
				vm(j)=temp;
%
				itemp=vrank(i);
                vrank(i)=vrank(j);  
				vrank(j)=itemp;
            end    
			if(dm(i)<dm(j))
				temp=dm(i);
                dm(i)=dm(j);  
				dm(j)=temp;
%
			    itemp=drank(i);
                drank(i)=drank(j);  
				drank(j)=itemp;
            end    
			if(em(i)<em(j))
				temp=em(i);
                em(i)=em(j);  
				em(j)=temp;
%
			    itemp=erank(i);
                erank(i)=erank(j);  
				erank(j)=itemp;
            end
			if(im(i)<im(j))
				temp=im(i);
                im(i)=im(j);  
				im(j)=temp;
%
			    itemp=irank(i);
                irank(i)=irank(j);  
				irank(j)=itemp;
            end
			if(cm(i)<cm(j))
				temp=cm(i);
                cm(i)=cm(j);  
				cm(j)=temp;
%
			    itemp=crank(i);
                crank(i)=crank(j);  
				crank(j)=itemp;
            end               
			if(km(i)<km(j))
				temp=km(i);
                km(i)=km(j);  
				km(j)=temp;
%
			    itemp=krank(i);
                krank(i)=krank(j);  
				krank(j)=itemp;
            end   
 			if(dskm(i)<dskm(j))
				temp=dskm(i);
                dskm(i)=dskm(j);  
				dskm(j)=temp;
%
			    itemp=dskrank(i);
                dskrank(i)=dskrank(j);  
				dskrank(j)=itemp;
            end              
        end
    end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    for i=1:rntrials
        for j=1:rntrials
            if(yrank(i)==j)
				pyrank(j)=i;
                break;
            end   
        end
    end
%
    for i=1:rntrials
        for j=1:rntrials
            if(vrank(i)==j)
				pvrank(j)=i;
                break;
            end   
        end
    end
%
    for i=1:rntrials
        for j=1:rntrials
            if(drank(i)==j)
				pdrank(j)=i;
                break;
            end
        end
    end
%
    for i=1:rntrials
        for j=1:rntrials
            if(erank(i)==j)
				perank(j)=i;
                break;
            end   
        end
    end
%
    for i=1:rntrials
        for j=1:rntrials
            if(irank(i)==j)
				pirank(j)=i;
                break;
            end   
        end
    end
%
    for i=1:rntrials
        for j=1:rntrials
            if(crank(i)==j)
				pcrank(j)=i;
                break;
            end   
        end
    end    
%
    for i=1:rntrials
        for j=1:rntrials
            if(krank(i)==j)
				pkrank(j)=i;
                break;
            end   
        end
    end
%
    for i=1:rntrials
        for j=1:rntrials
            if(dskrank(i)==j)
				pdskrank(j)=i;
                break;
            end   
        end
    end
%
    nmax=0.;
    iwin=0;
%


%    out1=sprintf('\n iw=%d ew=%d dw=%d vw=%d aw=%d  \n',iw,ew,dw,vw,aw);
%    disp(out1);
    
%    pirank
    
    try
        nrank=((iw*pirank+ew*perank)+(dw*pdrank+vw*pvrank+aw*pyrank)+cw*pcrank+kw*pkrank+dskw*pdskrank);
    catch
        warndlg('nrank failed');
        return;
    end
        
%    for i=1:rntrials 
%        out1=sprintf('i=%d  nrank=%d  pyrank=%d',i,nrank(i),pyrank(i));
%        disp(out1);
%    end
 
    
%
    for i=1:rntrials 
%
		if( nrank(i)>nmax)
%%            out1=sprintf('p1 %d  %8.4g  %8.4g ',i,dm(i),displacement_limit);
%%            disp(out1);
			nmax=nrank(i);
			iwin=i;
        end
    end
%
    nmax=0.;
    for i=1:rntrials 
%
		if( nrank(i)>nmax && dm(pdrank(i))<=displacement_limit)
%%         out1=sprintf('p2 %d  %8.4g  %8.4g ',i,dm(i),displacement_limit);
%%        disp(out1);
			nmax=nrank(i);
			iwin=i;
        end
    end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%
