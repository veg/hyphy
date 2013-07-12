class Mock_AVLList : public _AVLList {
 public:
  MOCK_METHOD1(Clear,
      void());
  MOCK_METHOD1(HasData,
      bool(long));
  MOCK_METHOD1(ReorderList,
      void());
  MOCK_METHOD3(InsertData,
      long(BaseRef, long, bool));
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD3(Traverser,
      long(, , ));
  MOCK_METHOD0(GetRoot,
      long(void));
  MOCK_METHOD1(DeleteXtra,
      void(long));
  MOCK_METHOD1(DeleteAll,
      void(bool cL));
};
