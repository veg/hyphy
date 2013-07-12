class Mock_List : public _List {
 public:
  MOCK_METHOD1(GetItem,
      BaseRef(const unsigned long));
  MOCK_METHOD1(=,
      operator _List(_List));
  MOCK_METHOD1(BinaryFind,
      long(BaseRef));
  MOCK_METHOD2(Compare,
      long(long, long));
  MOCK_METHOD2(Compare,
      long(BaseRef, long));
  MOCK_METHOD1(Clear,
      void());
  MOCK_METHOD1(Duplicate,
      void(const BaseRef));
  MOCK_METHOD1(DeleteList,
      void(const _SimpleList));
  MOCK_METHOD2(Find,
      long(, long));
  MOCK_METHOD2(FindPointer,
      long(BaseRef, long));
  MOCK_METHOD4(FindString,
      long(, long, bool, long));
  MOCK_METHOD1(FreeUpMemory,
      long(long));
  MOCK_METHOD3(InsertElement,
      void(BaseRef, long, bool));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD1(toFileStr,
      void(FILE));
};
