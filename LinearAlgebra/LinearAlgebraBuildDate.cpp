const char* LinearAlgebraBuildDate();

const char* LinearAlgebraBuildDate()
{
#ifdef BUILD_DATE
   return BUILD_DATE;
#else
   return "Unknown";
#endif
}

