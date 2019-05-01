//
//  l2_enet.cpp
//  osc
//
//  Created by Futao Zhang on 5/11/2018.
//  Copyright © 2018 Futao Zhang. All rights reserved.
//

#include "l2_enet.hpp"
namespace ELNET
{
    void  get_int_parms(double &sml,double &eps,double &big,int &mnlam,double &rsqmax,double &pmin,double &exmx)
    {
        double sml0=1e-5, eps0=1e-6, big0=9.9e35, rsqmax0=0.999,pmin0=1e-9, exmx0=250.0;
        int mnlam0=5;
        sml=sml0;
        eps=eps0;
        big=big0;
        mnlam=mnlam0;
        rsqmax=rsqmax0;
        pmin=pmin0;
        exmx=exmx0;
    }
    void elnet1(double &parm,int &ni,VectorXi &ju,VectorXd &vp,MatrixXd &cl,VectorXd &g,int &no,int &ne,int &nx,MatrixXd &x,int &nlam, double &flmin,VectorXd &ulam,double &thr,int &maxit,VectorXd &xv,int &lmu,MatrixXd &ao,VectorXi &ia,VectorXi &kin,VectorXd &rsqo,VectorXd &almo,int &nlp,int &jerr)
    {/*
        MatrixXd c(ni, nx);
        double sml,eps,big,rsqmax,pmin,exmx;
        int mnlam;
        get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx);
        VectorXd a=VectorXd::Zero(ni), da(ni);
        VectorXi mm=VectorXi::Zero(ni);
        double bta=parm;
        double omb=1-bta;
        double alf=1.0;
        if(flmin<1.0)
        {
            double eqs=max(eps,flmin);
            alf=pow(eqs,1.0/(nlam-1));
        }
        double rsq=0.0;
        nlp=0;
        double nin=nlp;
        double iz=0;
        double mnl= min(mnlam,nlam);
        double alm=0;
        for(int m=0;m<nlam;m++)
        {
            if(flmin<1.0)
            {
                if(m<=2)
                {
                    if(m!=1)
                    {
                        alm=0.0;
                        for(int j=0;j<ni;j++)
                        {
                            if(ju[j]!=0 && vp[j]>0)
                            {
                                alm=max(alm,abs(g(j))/vp(j));
                            }
                            alm=alf*alm/max(bta,1.0e-3);
                        }
                    }
                    else
                    {
                        alm=big;
                    }
                }
                else
                {
                    alm=alm*alf;
                }
            }
            else
            {
                alm=ulam(m);
            }
            double dem=alm*omb;
            double ab=alm*bta;
            double rsq0=rsq;
            double jz=1;
            while()
            {
                double dlx=0.0;
                if(iz*jz==0)
                {
                    nlp=nlp+1;
                    
                    for(int k=0;k<ni;k++)
                    {
                        if(ju[k]!=0)
                        {
                            double ak=a(k);
                            double u=g(k)+ak*xv(k);
                            double v=abs(u)-vp(k)*ab;
                            a(k)=0.0;
                            if(v>0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)));
                            if(a(k)!=ak && nx<nin)
                            {
                                if(mm(k)==0)
                                {
                                    nin=nin+1;
                                        for(int j=0;j<ni;j++)
                                        {
                                            if(ju[j]!=0)
                                            {
                                                if(mm(j) == 0)
                                                {
                                                    if(j != k)
                                                    {
                                                        c(j,nin)=x.col(j).dot(x.col(k));
                                                    }
                                                    else
                                                    {
                                                        c(j,nin)=xv(j);
                                                    }
                                                }
                                                else
                                                {
                                                    c(j,nin)=c(k,mm(j));
                                                }
                                            }
                                        }
                                        mm(k)=nin;
                                        ia(nin)=k;
                                    
                                }
                                double del=a(k)-ak ;
                                double rsq=rsq+del*(2.0*g(k)-del*xv(k));
                                dlx=max(xv(k)*del*del,dlx);
                                for(int j=0;j<ni;j++)
                                {
                                    if(ju(j)!=0) g(j)=g(j)-c(j,mm(k))*del;
                                }
                            }
                        }
                    }
                    if(dlx < thr) goto 10352;
                    if(nin > nx) goto 10352;
                    if(nlp > maxit){
                        jerr=-m;
                        return;
                    }
                }
                iz=1;
                for(int j=0;j<nin;j++) da[j]=a(ia[j]);
                while(dlx>=thr)
                {
                    nlp=nlp+1;
                    dlx=0.0;
                    for(int l=0;l<nin;l++)
                    {
                        int k=ia(l);
                        double ak=a(k);
                        double u=g(k)+ak*xv(k);
                        double v=abs(u)-vp(k)*ab;
                        a(k)=0.0;
                        if(v>0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)));
                        if(a(k)!=ak)
                        {
                            double del=a(k)-ak;
                            rsq=rsq+del*(2.0*g(k)-del*xv(k));
                            dlx=max(xv(k)*del*del,dlx);
                            for(int j=1;j<nin;j++)
                            {
                                g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del;
                                continue
                            }
                        }
                    }
                    if(nlp > maxit){
                        jerr=-m;
                        return;
                    }
                }
                for(int j=0;j<nin;j++)  da[j]=a(ia[j])-da[j];
                for(int j=0;j<ni;j++) {
                    if(mm(j)==0){
                        if(ju(j)!=0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin));
                    }
                }
                
                jz=0;
            }
           
            
        }
      */
    }
    void standard(int &no,int &ni,MatrixXd &x,VectorXd &y,VectorXd &w,int &isd,int &intr,VectorXi &ju,VectorXd &g,VectorXd &xm,VectorXd &xs,double &ym,double &ys,VectorXd &xv,int &jerr)
    {
        VectorXd v(no);
        double sumw=0.0;
        for(int i=0;i<no;i++) sumw+=w[i];
        for(int i=0;i<no;i++)
        {
            w[i]/=sumw;
            v[i]=sqrt(w[i]);
        }
        if(intr) //intercept
        {
            for(int i=0;i<ni;i++)
            {
                if(ju[i]!=0)
                {
                    xm(i)=w.dot(x.col(i));
                    x.col(i)=x.col(i).array()-xm(i);
                    x.col(i)=v.array()*x.col(i).array();
                    xv(i)=x.col(i).dot(x.col(i));
                    if(isd) xs(i)=sqrt(xv(i));
                }
            }
            if(isd)
            {
                for(int i=0;i<ni;i++)
                {
                    if(ju[i]!=0)
                    {
                        x.col(i)=x.col(i)/xs(i);
                    }
                }
                for(int i=0;i<ni;i++) xv[i]=1.0;
            }
            else
            {
                for(int i=0;i<ni;i++) xs[i]=1.0;
            }
            ym=w.dot(y);
            y=y.array()-ym;
            y=v.array()*y.array();
            ys=sqrt(y.dot(y));
            y=y/ys;
        }
        else
        {
            ym=0;
            y=v.array()*y.array();
            double vy=v.dot(y);
            ys=sqrt(y.dot(y)-vy*vy);
            y/=ys;
            for(int i=0;i<ni;i++)
            {
                if(ju[i]!=0)
                {
                    xm(i)=0;
                    x.col(i)=v.array()*x.col(i).array();
                    xv(i)=x.col(i).dot(x.col(i));
                    if(isd)
                    {
                        double xbq=v.dot(x.col(i));
                        xbq*=xbq;
                        double vc=xv(i)-xbq;
                        xs(i)=sqrt(vc);
                        x.col(i)=x.col(i)/xs(i);
                        xv(i)=1+xbq/vc;
                    } else {
                        xs(i)=1.0;
                    }
                }
            }
        }
        for(int i=0;i<ni;i++) g[i]=0.0;
        for(int i=0;i<ni;i++)
        {
            if(ju[i]!=0) g[i]=y.dot(x.col(i));
        }
    }
    void chkvars(int &no, int &ni,MatrixXd &x,VectorXi &ju)
    {
        for(int j=0;j<ni;j++)
        {
            ju[j]=0;
            double t=x(0,j);
            for(int i=1;i<no;i++)
                if(x(i,j)!=t) ju[j]=1;
        }
    }
    void elnetu(double &parm,int &no, int &ni,MatrixXd &x, VectorXd &y,VectorXd &weights,VectorXi &jd,VectorXd &vp,MatrixXd &cl,int &ne,int &nx,int &nlam,double &flmin,VectorXd &ulam,double &thr,int &isd,int &intr,int &maxit,  int &lmu,VectorXd &a0,MatrixXd &ca,VectorXi &ia,VectorXi &nin,VectorXd &rsq,VectorXd &alm,int &nlp,int &jerr )
    {
        VectorXd xm(ni), xs(ni),g(ni),xv(ni),vlam(nlam);
        VectorXi ju(ni);
        double ym, ys;
        chkvars(no, ni, x, ju);
        if(jd[0]>0)
            for(int i=1;i<jd[0]+1;i++) ju[jd[i]]=0;
        if(ju.maxCoeff()<=0)
        {
            jerr=7777;
            return;
        }
        standard(no,ni,x,y,weights,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr);
        cl=cl/ys;
        if(isd>0)
        {
            for(int i=0;i<ni;i++)
                cl.col(i)=cl.col(i).array()*xs(i);
        }
        if(flmin>=1.0) vlam=ulam/ys;
        elnet1(parm,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr);
        
    }
    void elnet(MatrixXd &x, vector<double> &y,vector<double> &weights, vector<double> &offset, char* typeGau, double alpha, int nobs, int nvars, vector<int> &jd,vector<double> &vp,MatrixXd &cl,int ne,int nx,int nlam,double flmin,vector<double> &ulam,double thresh,int isd,int intr,char* vnames,int maxit)
    {
        //x(no,ni),y(no),w(no),vp(ni),cl(2,ni),ulam(nlam)
        
        int ka=1; // 1 for covariance 2 for naive
        if(offset.size()>0)
            for(int i=0;i<y.size();i++)
                y[i]=y[i]-offset[i];
        double ybar=weight_mean(y, weights);
        double nulldev=0.0;
        for(int i=0;i<y.size();i++)
            nulldev+=weights[i]*(y[i]-ybar)*(y[i]-ybar);
        if(nulldev==0)
        {
            LOGPRINTF("ERROR: y is constant.\n");
            TERMINATE();
        }
        double parm=alpha;
        int no=nobs;
        int ni=nvars;
        int lmu=1;
        vector<double> a0(nlam);
        MatrixXd ca(nx,nlam);
        vector<int> ia(nx);
        vector<int> nin(nlam);
        vector<double> rsq(nlam);
        vector<double> alm(nlam);
        int nlp=0;
        int jerr=0;
        vector<double> vq(ni);
        if(*max_element(vp.begin(), vp.end())<=0.0)
        {
            jerr=10000;
            return;
        }
        double sumvq=0.0;
        for(int i=0;i<ni;i++)
        {
            vq[i]=vp[i]>0?vp[i]:0;
            sumvq+=vq[i];
        }
        sumvq=ni/sumvq;
        for(int i=0;i<ni;i++) vq[i]*=sumvq;
       // if(ka==1) elentu();
        //else eletn();
        
    }
    void glmnet(MatrixXd &x,VectorXd &y,char* family,double alpha,int nlam,bool standardize, bool intercept)
    {
        //alpha=1, nlam=100 as default
        long nobs=x.rows(), nvars=x.cols(),nrowy=y.size();
        if(nobs!=nrowy)
        {
            LOGPRINTF("number of observations in y (%ld) not equal to the number of rows of x (%ld)",nobs,x.rows());
            TERMINATE();
        }
        VectorXd weight=VectorXd::Ones(nobs), penaltyFacotr=VectorXd::Ones(nobs);
        long dfmax=nvars+1, pmax=min(dfmax*2+20,nvars);
        /*
        lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4)
        lambda=NULL;standardize=TRUE;intercept=TRUE;thresh=1e-7;dfmax=nvars+1;pmax=min(dfmax*2+20,nvars);exclude=NULL;penalty.factor=rep(1,nvars);lower.limits=-Inf;upper.limits=Inf;maxit=100000;type.gaussian=ifelse(nvars<500,"covariance","naive");
        type.logistic=c("Newton","modified.Newton");standardize.response=FALSE;type.multinomial=c("ungrouped","grouped")
        //glmnet.control()
        fdev=1.0e-5, devmax=0.999, eps=1.0e-6, big=9.9e35, mnlam=5, pmin=1.0e-9,
        exmx=250.0,prec=1e-10,mxit=100,factory=FALSE
        is.sparse=FALSE
        //jd: jd[0] is the total number of variable to be excluded. the elements onwards are the index of the variable to be excluded. if exclude=NULL, then jd=0
        */
    }
    void cv_glmnet(MatrixXd &x,VectorXd &y, char* measure, VectorXd &lambda,VectorXd &dev, int nfolds, VectorXd &cvm, VectorXd &cvsd,VectorXi &nzero,double &lambda_min, double &lambda_1se)
    {
        int weigth=1;
        char* offset=NULL;
        bool grouped=true;
        //glmnet();
        
    }
}
