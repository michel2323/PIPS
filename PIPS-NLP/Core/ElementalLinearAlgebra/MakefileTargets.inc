all_stoch_libs: $(libooqpemtl_STATIC)
all_stoch_headers: stoch_headers

stoch_headers: $(EMTLHEADERS)

$(libooqpemtl_STATIC) : $(libooqpemtl_OBJ)
	ar ru $@ $(libooqpemtl_OBJ)
	$(RANLIB) $@

$(localincludedir)/%.h : $(EMTLLINALGDIR)/%.h
	$(INSTALL) -m 444 $< $@
