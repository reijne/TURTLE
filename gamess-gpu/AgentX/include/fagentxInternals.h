
fortint axcwrapproxy(fortintc* proxylen, char* proxy);
fortint axcwrapgeturi(fortintc* urilen, char* uri);
fortint axcwrapdatageturi(fortintc* urilen, char* uri);
fortint axcwrapselect(fortintc* entitylen, char* entity);
fortint axcwrapselectnext();
fortint axcwrapdeselect();
fortint axcwrapbaseuri(fortintc* local_base_uri_stringlen, char* local_base_uri_string);
fortint axcwrapparserfinish();
fortint axcwrapparserstart();
fortint axcwraprefine(fortintc* expressionlen, const char* expression);
fortint axcwrapformat(fortint *dataFormat);
fortint axcwrapcurrent();
fortint axcwraplen();
void axcwrapvalue(fortintc* strvarlen, char* strvar);
fortint axcwrapbuffer(fortint* size);
fortint axcwrapnamespace(fortintc* namespacefURIlen, char* namespacefURI);
fortint axcwrapcache(fortint *cachesize);
fortint axcwrapselectno(fortint *no);
