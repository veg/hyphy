class Mock_Trie : public _Trie {
 public:
  MOCK_METHOD0(toStr,
      BaseRef(void));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef storage));
  MOCK_METHOD1(Clear,
      void(bool));
};
