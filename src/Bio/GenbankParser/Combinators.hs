module Bio.GenbankParser.Combinators (
    eol
  , whiteSpace
  , notChar
  , noEol
  , noWhiteSpace
)	where

import Text.ParserCombinators.Parsec
import Text.ParserCombinators.Parsec.Token hiding (whiteSpace)
import Text.ParserCombinators.Parsec.Language (emptyDef)    
import Control.Applicative hiding ((<|>),(<?>))

eol = try (string "\n")
    <|> try (string "\r")
    <|> try (string "\n\r")
    <?> "End of Line"

whiteSpace = string " "
            <|> string "\t"
            <|> eol
            <?> "Whitespace"

notChar x = notFollowedBy x *> anyChar

noWhiteSpace = notChar whiteSpace

noEol = notChar eol
