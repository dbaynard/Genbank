{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell #-}
module Test (
)	where

import Text.Parsec.Prim (runP)
import Data.Text.Lazy ()
{-import qualified Data.Text.Lazy as DT (unwords)-}
import qualified Data.ByteString.Lazy.Char8 as B (pack)
import qualified Data.Text.Lazy.IO as DTI
import Text.Parsec hiding (space, spaces)
import qualified Text.Parsec.Text.Lazy as T
import qualified Text.Parsec.ByteString.Lazy as P ()
import Control.Applicative hiding ((<|>), many, optional)
{-import Text.Parsec.Indent-}
import Data.Monoid
import Control.Monad hiding (sequence)
import Data.Traversable (sequence, sequenceA)

import Control.Lens hiding (locus)

import Bio.GenbankData

import Language.Haskell.TH
makeLensesWith ?? ''Genbank $ lensRules
    & lensField .~ \_ _ name -> [TopName (mkName $ nameBase name ++ "L")]

-- | @parseFromFile p filePath@ runs a lazy bytestring parser @p@ on the
-- input read from @filePath@ using 'ByteString.Lazy.Char8.readFile'. Returns either a 'ParseError'
-- ('Left') or a value of type @a@ ('Right').
--
-- >  main    = do{ result <- parseFromFile numbers "digits.txt"
-- >              ; case result of
-- >                  Left err  -> print err
-- >                  Right xs  -> print (sum xs)
-- >              }
parseFromFile :: T.Parser a -> String -> IO (Either ParseError a)
parseFromFile p fname
    = do input <- DTI.readFile fname
         return (runP p () fname input)

{-type IParser a = IndentParser Text () a-}
{-
 -parseFromFile' :: IParser a -> String -> IO (Either ParseError a)
 -parseFromFile' p fname
 -    = do input <- DTI.readFile fname
 -         return (runParserT p () fname input)
 -}

(~-~) :: (Applicative f) => ASetter s t a b -> f b -> s -> f t
(~-~) = f
    where
        lens `f` parser = ( (lens .~) <$> parser <*> ) . pure
infixr 5 ~-~

space :: T.GenParser st Char
space = oneOf "         "
        <?> "Space"

spaces = many space
        <?> "Zero or more spaces"

spaces1 = many1 space
        <?> "One or more space"

smallIndent = space *> space

middlIndent = smallIndent *> smallIndent *> space

largeIndent = middlIndent *> many1 space

whiteSpace = space
            <|> char '\t'
            <|> endOfLine
            <?> "Whitespace"

notChar x = notFollowedBy x *> anyChar

noWhiteSpace = notChar whiteSpace

noEol = notChar endOfLine

itemJoin = try $ endOfLine *> lookAhead noWhiteSpace

subItemJoin = try $ endOfLine *> smallIndent

featureJoin = try $ endOfLine *> middlIndent

gLocusParser x = do
        string "LOCUS"
        spaces1
        locus <- many1 noWhiteSpace
        spaces1
        length <- many1 noWhiteSpace 
        string " bp"
        spaces1
        moleculeType <- many1 noWhiteSpace 
        spaces1
        circular <- optionMaybe . try . choice . fmap string $ ["linear", "circular", "LINEAR", "CIRCULAR"]
        spaces
        division <- optionMaybe . choice . fmap (try . string) $ ["PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT", "VRL", "PHG", "SYN", "UNA", "EST", "PAT", "STS", "GSS", "HTG", "HTC", "ENV"]
        spaces
        creationDate <- many1 noEol
        itemJoin
        {-return (locus, length, moleculeType, circular, division, creationDate)-}
        return $ x { locus = B.pack locus
                   , genbankLength = read length
                   , moleculeType = B.pack moleculeType
                   , circular = B.pack <$> circular
                   , division = B.pack <$> division
                   , creationDate = B.pack creationDate
                   }

gDefinitionParser = gSimpleAttribParser "DEFINITION"

gDefinitionParser' modifier  = setter modifier `fmap` gSimpleAttribParser "DEFINITION"

setter modifier item record = record { definition = modifier B.pack item }
{-
 -gDefinitionParser' x = setter definition id x `fmap` gSimpleAttribParser "DEFINITION"
 -
 -setter item modifier record = record { definition = modifier B.pack item }
 -}

-- gLocusParser defGenbank <**> fmap (\y x -> x { definition = B.pack y } ) gDefinitionParser

gAccessionParser = gSimpleAttribParser "ACCESSION"

gVersionParser = gSimpleAttribParser "VERSION"

gKeywordsParser = gSimpleAttribParser "KEYWORDS"

gCommentParser = gSimpleAttribParser "COMMENT"

gAttribParsers = [ gDefinitionParser
                 , gAccessionParser
                 , gVersionParser
                 , gKeywordsParser
                 , gCommentParser
                 , gSourceParser
                 , gReferenceParser
                 ]

gFundamentalAttribParser end identifier = string identifier <* spaces1
        *> multiLine
        <* end

gSimpleAttribParser = gFundamentalAttribParser itemJoin

gSubAttribParser = gFundamentalAttribParser subItemJoin

gFeatureAttribParser = gFundamentalAttribParser featureJoin

gAmbigAttribParser = gFundamentalAttribParser (subItemJoin <|> itemJoin)

multiLineGen combinator separator = combinator
        <$> many1 noEol `sepBy` try separator

multiLine = multiLineGen unwords (endOfLine *> largeIndent)

gSourceParser = f
        <$> gSubAttribParser "SOURCE" <*> gOrganismParser
        where
            f = (<>) . (<> " ")
            gOrganismParser = gSimpleAttribParser "ORGANISM"

gReferenceParser = f
        <$> gAmbigAttribParser "REFERENCE"
        <*> h 
        where
            f = (<>) . (<> " ")
            h = unwords <$> g gReferenceAttribs
            g = many . try . choice . fmap gAmbigAttribParser
            gReferenceAttribs = [ "AUTHORS"
                                , "TITLE"
                                , "JOURNAL"
                                , "PUBMED" ]

gFeatureTableParser = gFeatureAttribParser "FEATURES"
                    *> gFeatureParser `sepBy` try featureJoin
                    <* itemJoin

gFeatureParser = (,,)
                <$> many1 noWhiteSpace <* spaces1
                <*> gSpanParser <* featureSep
                <*> f `sepBy` try featureSep
                where
                    f = (,) <$> manyTill noWhiteSpace (char '=') <*> multiLineFeat
                    g = endOfLine *> largeIndent
                    featureSep = g *> char '/'
                    featureIntSep = g *> notFollowedBy (char '/')
                    multiLineFeat = multiLineGen (foldr (<>) mempty) featureIntSep

-- TODO
gSpanParser = many1 noEol

{-
 -gSpanParser = try f
 -            <|> g
 -            <?> "coordinates"
 -        where
 -            f = 
 -            g =
 -                <$> optionMaybe (char '<') *> many1 digit
 -                <*> string ".." *> many1 digit <* optionMaybe (char '>')
 -            parentheses = between (char '(') (char ')')
 -}

gOriginParser = join <$> nucline
        where
            nucline = string "ORIGIN" *> endOfLine *> many1 nucleotides <* string "//"
            nucleotides = join <$> nucString'
            nucString' = spaces *> many1 digit *> many1 (try nucString) <* spaces <* endOfLine
                  <?> "numbered line of nucleotide strings"
            nucString = spaces *> (many1 . oneOf $ nucLetters)
                  <?> "space delimited nucleotide sequence"
            nucLetters = "acgnturykmswbdhvACGNTURYKMSWBDHV"

{-(\x -> parseFromFile x "sample.gb") -}
testParse = (,,,) <$> gLocusParser defGenbank <*> many (choice gAttribParsers) <*> gFeatureTableParser <*> gOriginParser
