(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (or A13 A3))) (let ((.def_1 (not .def_0))) (let ((.def_2 (= A0 A11))) (let ((.def_3 (and .def_2 .def_1))) (let ((.def_4 (not .def_3))) (let ((.def_5 (or A7 A3))) (let ((.def_6 (not A19))) (let ((.def_7 (or A22 .def_6))) (let ((.def_8 (= .def_7 .def_5))) (let ((.def_9 (and .def_8 .def_4))) (let ((.def_10 (not .def_9))) (let ((.def_11 (not A3))) (let ((.def_12 (or A6 .def_11))) (let ((.def_13 (not A21))) (let ((.def_14 (and .def_13 A16))) (let ((.def_15 (and .def_14 .def_12))) (let ((.def_16 (not .def_15))) (let ((.def_17 (not A13))) (let ((.def_18 (= A24 .def_17))) (let ((.def_19 (not .def_18))) (let ((.def_20 (and A22 A10))) (let ((.def_21 (and .def_20 .def_19))) (let ((.def_22 (= .def_21 .def_16))) (let ((.def_23 (and .def_22 .def_10))) (let ((.def_24 (not .def_23))) (let ((.def_25 (and .def_17 A18))) (let ((.def_26 (and A7 A11))) (let ((.def_27 (not .def_26))) (let ((.def_28 (and .def_27 .def_25))) (let ((.def_29 (not A18))) (let ((.def_30 (= A9 .def_29))) (let ((.def_31 (not A17))) (let ((.def_32 (and .def_31 A15))) (let ((.def_33 (not .def_32))) (let ((.def_34 (and .def_33 .def_30))) (let ((.def_35 (not .def_34))) (let ((.def_36 (or .def_35 .def_28))) (let ((.def_37 (not A12))) (let ((.def_38 (or A16 .def_37))) (let ((.def_39 (or A2 A15))) (let ((.def_40 (and .def_39 .def_38))) (let ((.def_41 (not A16))) (let ((.def_42 (or A17 .def_41))) (let ((.def_43 (not A6))) (let ((.def_44 (or A5 .def_43))) (let ((.def_45 (or .def_44 .def_42))) (let ((.def_46 (not .def_45))) (let ((.def_47 (and .def_46 .def_40))) (let ((.def_48 (and .def_47 .def_36))) (let ((.def_49 (= .def_48 .def_24))) (let ((.def_50 (not .def_49))) (let ((.def_51 (and A4 A8))) (let ((.def_52 (not .def_51))) (let ((.def_53 (= .def_32 .def_52))) (let ((.def_54 (not .def_53))) (let ((.def_55 (not A7))) (let ((.def_56 (and .def_55 .def_43))) (let ((.def_57 (not .def_56))) (let ((.def_58 (not A11))) (let ((.def_59 (not A5))) (let ((.def_60 (and .def_59 .def_58))) (let ((.def_61 (not .def_60))) (let ((.def_62 (and .def_61 .def_57))) (let ((.def_63 (or .def_62 .def_54))) (let ((.def_64 (not .def_63))) (let ((.def_65 (= .def_13 A6))) (let ((.def_66 (not .def_65))) (let ((.def_67 (and .def_11 A12))) (let ((.def_68 (not .def_67))) (let ((.def_69 (and .def_68 .def_66))) (let ((.def_70 (= A12 A1))) (let ((.def_71 (not .def_70))) (let ((.def_72 (and .def_6 .def_31))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and .def_73 .def_71))) (let ((.def_75 (or .def_74 .def_69))) (let ((.def_76 (and .def_75 .def_64))) (let ((.def_77 (not .def_76))) (let ((.def_78 (not A0))) (let ((.def_79 (or .def_78 A7))) (let ((.def_80 (not A15))) (let ((.def_81 (and A14 .def_80))) (let ((.def_82 (not .def_81))) (let ((.def_83 (and .def_82 .def_79))) (let ((.def_84 (not .def_83))) (let ((.def_85 (= .def_43 A24))) (let ((.def_86 (not .def_85))) (let ((.def_87 (and A4 A15))) (let ((.def_88 (and .def_87 .def_86))) (let ((.def_89 (and .def_88 .def_84))) (let ((.def_90 (not A4))) (let ((.def_91 (or .def_90 .def_58))) (let ((.def_92 (not A22))) (let ((.def_93 (or .def_92 .def_13))) (let ((.def_94 (or .def_93 .def_91))) (let ((.def_95 (not A2))) (let ((.def_96 (and A15 .def_95))) (let ((.def_97 (or A7 A21))) (let ((.def_98 (not .def_97))) (let ((.def_99 (or .def_98 .def_96))) (let ((.def_100 (and .def_99 .def_94))) (let ((.def_101 (or .def_100 .def_89))) (let ((.def_102 (and .def_101 .def_77))) (let ((.def_103 (not .def_102))) (let ((.def_104 (or .def_103 .def_50))) (let ((.def_105 (not .def_104))) .def_105)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)