def scanpy_workflow(ad,pathway,batch_key='sampleID',span=0.3):
    import scanpy as sc
    ad=sc.read(ad)
    ad.layers['counts']=ad.X.copy()
    sc.pp.normalize_total(ad,target_sum=1e4)
    sc.pp.log1p(ad)
    ad.raw=ad
    sc.pp.highly_variable_genes(
        ad,
        flavor='seurat_v3',
        layer='counts',
        n_top_genes=2000,
        batch_key=batch_key,
        subset=False,
        span=span
    
        
        
    )
    sc.tl.pca(ad, svd_solver='arpack', use_highly_variable=True)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad)
    ad.obsm['X_umapRaw'] = ad.obsm['X_umap']
    sc.external.pp.harmony_integrate(ad, key='sampleID', basis='X_pca', adjusted_basis='X_pca_harmony',max_iter_harmony=20)

    y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF_assessment(ad.obsm['X_pca_harmony'], ad.obs['cellType'],n=500)
    prob=pd.DataFrame(y_prob,columns=clf.classes_)
    pred=pd.DataFrame(y_pred,columns=['predict'])
    sccaf=pd.concat([pred,prob],axis=1)
    sccaf.to_csv(pathway+'/'+'sccaf-result.csv')
    y_test.to_csv(pathway+'/'+'y_test-result.csv')
